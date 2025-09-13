### A Pluto.jl notebook ###
# v0.19.46

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

# ╔═╡ 1cfa1e53-7cf7-4fa2-a36a-e9d36b134cad
bin = -5000:50:5000

# ╔═╡ dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
begin

	pplt.close(); fig,axs = pplt.subplots(nrows=2,ncols=6,sharey=0,axwidth=1)
	
	for ax in axs
		ax.format(xlabel=L"$<h'>$ / kJ hPa",xlim=(-5000,5000),ylim=(-5000,5000))
	end
	axs[1].format(ylabel=L"$<c_pT'>$ / kJ hPa",ultitle=L"$a_m = 0.1$ day$^{-1}$")
	axs[2].format(ultitle=L"$a_m = 0.2$ day$^{-1}$")
	axs[3].format(ultitle=L"$a_m = 0.5$ day$^{-1}$")
	axs[4].format(ultitle=L"$a_m = 1$ day$^{-1}$")
	axs[5].format(ultitle=L"$a_m = 2$ day$^{-1}$")
	axs[6].format(ultitle=L"$a_m = 5$ day$^{-1}$")
	axs[7].format(ylabel=L"$<L_vq'>$ / kJ hPa")

	for ii in vcat(2:6,8:12)
		axs[ii].format(yformatter="none")
	end

	for ic = 1 : 6

		if wtgscheme == "DGW"
			config = dampingstrprnt(configWTG[ic+6])
		else
			config = relaxscalestrprnt(configWTG[ic+8])
		end

		ds = NCDataset(datadir("moisturemode","$wtgscheme-$expname-$config.nc"))
		hp = (ds["hp"][:,11:15])
		Tp = (ds["Tp"][:,11:15]) * 1.0035
		qp = (ds["qp"][:,11:15])/1000 * 2500.9
		close(ds)

		hp = (hp .- mean(hp[1:100,:],dims=1))[:]
		Tp = (Tp .- mean(Tp[1:100,:],dims=1))[:]
		qp = (qp .- mean(qp[1:100,:],dims=1))[:]
	
		axs[ic].scatter(hp,Tp,s=1,a=0.05)
		axs[ic+6].scatter(hp,qp,s=1,a=0.05)
		axs[ic].format(suptitle="$wtgscheme | $expname")

	end

	fig.savefig(plotsdir("09-moisturemode-$wtgscheme-$expname.png"),dpi=400)
	load(plotsdir("09-moisturemode-$wtgscheme-$expname.png"))
	
end

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
# ╠═a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╠═1cfa1e53-7cf7-4fa2-a36a-e9d36b134cad
# ╠═dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
# ╟─8fe5cd5d-28c5-4607-8dfb-83e4c0141a7c
