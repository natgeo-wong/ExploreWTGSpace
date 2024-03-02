### A Pluto.jl notebook ###
# v0.19.32

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
	using StatsBase
	using Trapz
	
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
	if wtgscheme == "DGW"
		configWTG = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	else
		configWTG = [
			0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
			1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
			10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
		]
	end
	nconWTG = length(configWTG)
	blues_WTG = pplt.get_colors("Blues",(nconWTG+2))
	lgd_WTG = Dict("frame"=>false,"ncols"=>3)
	md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
begin
	prcp = zeros(6000,15,nconWTG)
	gms  = zeros(6000,15,nconWTG)
	wdse = zeros(6000,15,nconWTG)
	wmse = zeros(6000,15,nconWTG)
	Γcrt = zeros(15,nconWTG)

	for ic in 1 : nconWTG

		if wtgscheme == "DGW"
			config = dampingstrprnt(configWTG[ic])
		else
			config = relaxscalestrprnt(configWTG[ic])
		end
		ds   = NCDataset(datadir("precipitation","DGW-T1282km300V64-$config.nc"))
		prcp[:,:,ic] = ds["precipitation"][:] / 24
		close(ds)
		ds  = NCDataset(datadir("gms","DGW-T1282km300V64-$config.nc"))
		gms[:,:,ic]  = ds["gms"][:]
		wdse[:,:,ic] = ds["wdse"][:]
		wmse[:,:,ic] = ds["wmse"][:]
		close(ds)

	end

	for ic in 1 : nconWTG, ien = 1 : 15
		Γcrt[ien,ic] = wdse[:,ien,ic] \ wmse[:,ien,ic]
	end

	rgms = sqrt.(wdse.^2 .+ wmse.^2)
	θgms = atan.(wmse,wdse) .- atan.(reshape(Γcrt,1,15,nconWTG))
	wdse_rot = rgms .* cos.(θgms)
	wmse_rot = rgms .* sin.(θgms)
	
	md"Loading precipitation and GMS data ..."
end

# ╔═╡ 9f2e48b9-ce3b-4803-931a-8e68a374a362
nt = 200

# ╔═╡ 0d882c1f-075d-4688-a1d6-0f5ec88b941f
smthtime = 12

# ╔═╡ 7f50f030-f2ec-4a45-a3df-5b1786a0d77a
md"Do Animation? $(@bind doanim PlutoUI.Slider(0:1,default=0))"

# ╔═╡ b194253b-a123-457e-800d-fbe44cea1adc
if isone(doanim)
	for it = 1 : (nt-(smthtime-1))
		pplt.close()
		fig,axs = pplt.subplots(axwidth=1.5,ncols=2,sharex=0,sharey=0)
	
		for icon = 1 : nconWTG
			axs[1].scatter(
				dropdims(mean(wdse_rot[6000-nt+it.+(0:(smthtime-1)),:,icon],dims=1),dims=1),
				dropdims(mean(wmse_rot[6000-nt+it.+(0:(smthtime-1)),:,icon],dims=1),dims=1),
				s=10,c=blues_WTG[icon+1]
			)
			axs[2].scatter(
				dropdims(mean(wmse_rot[6000-nt+it.+(0:(smthtime-1)),:,icon],dims=1),dims=1) ./
				dropdims(mean(wdse_rot[6000-nt+it.+(0:(smthtime-1)),:,icon],dims=1),dims=1),
				dropdims(mean(prcp[6000-nt+it.+(0:(smthtime-1)),:,icon],dims=1),dims=1) .-
				dropdims(mean(prcp[6000-nt-1+it.+(0:(smthtime-1)),:,icon],dims=1),dims=1),
				s=10,c=blues_WTG[icon+1]
			)
		end
	
		# axs[1].plot([-5,5],[-5,5]*0.4,c="k",linestyle="--",lw=0.5)
		# axs[1].text(0.9,0.32,L"$\Gamma_c \approx 0.4$",rotation=45,horizontalalignment="center", verticalalignment="bottom")
		axs[1].format(
			xlim=(-0.5,1.5),xlocator=-0.5:0.5:1.5,
			ylim=(-0.3,0.3),ylocator=-0.2:0.2:0.6,
			ylabel=L"<$w\cdot\partial_zh$>",
			xlabel=L"<$w\cdot\partial_zs$>",
		)
		
		axs[2].format(
			xlim=(-20,20),xscale="symlog",xlabel=L"GMS - $\Gamma_c$",
			ylim=(-0.2,0.2),yscale="symlog",ytickloc="r",ylabel=L"P / mm hr$^{-1}$",
			xscale_kw=Dict("linthresh"=>0.1),yscale_kw=Dict("linthresh"=>0.001),
		)
	
		mkpath(plotsdir("gmstest"))
		fig.savefig(plotsdir("gmstest","test-$it.png"),transparent=false,dpi=150)
	end
else
	md"Temporarily disable animation creation"
end

# ╔═╡ ee11960a-6124-4e17-8681-7e2b429e6248
load(plotsdir("gmstest","test-7.png"))

# ╔═╡ 8fe5cd5d-28c5-4607-8dfb-83e4c0141a7c
begin
	# pplt.close()
	# f1,a1 = pplt.subplots(axwidth=1.5,ncols=2,sharex=0,sharey=0)
	
	# for icon = 1 : nconWTG
	# 	a1[1].scatter(wdse[3601:end,:,icon],wmse[3601:end,:,icon],s=10,c=blues_WTG[icon+1])
	# 	a1[2].scatter(gms[3601:end,:,icon].-0.4,prcp[3601:end,:,icon],s=10,c=blues_WTG[icon+1])
	# end
	
	# a1[1].plot([-5,5],[-5,5]*0.4,c="k",linestyle="--",lw=0.5)
	# a1[1].text(0.9,0.32,L"$\Gamma_c \approx 0.4$",rotation=45,horizontalalignment="center", verticalalignment="bottom")
	# a1[1].format(
	# 	xlim=(-1,3),xlocator=-1:3,
	# 	ylim=(-0.4,1.2),ylocator=-0.4:0.4:1.2,
	# 	ylabel=L"<$w\cdot\partial_zh$>",
	# 	xlabel=L"<$w\cdot\partial_zs$>",
	# )
	
	# a1[2].format(
	# 	xlim=(-500,500),xscale="symlog",xlabel=L"GMS - $\Gamma_c$",
	# 	ylim=(0,2.5),yscale="symlog",ytickloc="r",ylabel=L"P / mm hr$^{-1}$",
	# 	yscale_kw=Dict("linthresh"=>0.1),
	# )
	
	# mkpath(plotsdir("gmstest"))
	# f1.savefig(plotsdir("gmstest","test.png"),transparent=false,dpi=150)
	# load(plotsdir("gmstest","test.png"))
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
# ╟─dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
# ╠═9f2e48b9-ce3b-4803-931a-8e68a374a362
# ╠═0d882c1f-075d-4688-a1d6-0f5ec88b941f
# ╟─7f50f030-f2ec-4a45-a3df-5b1786a0d77a
# ╟─b194253b-a123-457e-800d-fbe44cea1adc
# ╠═ee11960a-6124-4e17-8681-7e2b429e6248
# ╟─8fe5cd5d-28c5-4607-8dfb-83e4c0141a7c
