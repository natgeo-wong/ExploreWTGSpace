### A Pluto.jl notebook ###
# v0.20.19

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    #! format: off
    return quote
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
	using PlutoUI
	using NCDatasets, DSP
	using Printf, StatsBase
	using CairoMakie, MakieExtra, LaTeXStrings, ColorSchemes
	set_theme!(theme_latexfonts())
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
	"LGW" => "(LGW) Linear Gravity Waves [Kuang et al., 2008]",
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
	expname = "D0642km300V64"
md"**Experiment Set:** $expname"
end

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	if (wtgscheme == "DGW") || (wtgscheme == "LGW")
		configWTG = [0.02,0.05,0.1,0.5,1,2,5,10,20,50,100,200,500]
	else
		configWTG = [
			0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
			1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),
			10,10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
		]
	end
	nconWTG = length(configWTG)
	blues_WTG = get(ColorSchemes.Blues, range(0.0, 1.0, length=nconWTG))
	lgd_WTG = Dict("frame"=>false,"ncols"=>3)
	md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ dc06d1f6-b3ab-4c7e-9f30-cfe1bb5e9cbd
begin
	f1 = Figure()
	ax1 = Axis(
		f1[1,1],width=450,height=200,limits=(0,250,0,10),
		yscale=SymLog(0.01),yticks=[0,0.01,0.1,1,10],
		ylabel=L"Rainfall Rate / mm hr$^{-1}$",xlabel="Days"
	)

	for ic in 1 : nconWTG

		if (wtgscheme == "DGW") || (wtgscheme == "LGW")
			fnc = "$(dampingstrprnt(configWTG[ic])).nc"
		else
			fnc = "$(relaxscalestrprnt(configWTG[ic])).nc"
		end
		fnc = datadir("precipitation",wtgscheme,expname,fnc)
		if isfile(fnc)
			ds_dgwprcp = NCDataset(fnc)
			dt    = ds_dgwprcp["time"][:]
			prcp  = ds_dgwprcp["precipitation"][:,:] / 24
			for imem = 1 : 15
				idt = dropdims(mean(reshape(dt,24,:),dims=1),dims=1)
				iprcp = dropdims(mean(reshape(prcp[:,imem],24,:),dims=1),dims=1)
				lines!(ax1,idt,iprcp,color=blues_WTG[ic])
			end
			close(ds_dgwprcp)
		end

	end

	resize_to_layout!(f1)
	f1
end

# ╔═╡ 3f880613-b30e-417c-99e0-d78ea2302891
CairoMakie.save(plotsdir("02b-rce2wtg-$(expname)-$wtgscheme.png"),f1)

# ╔═╡ 0a74e728-1cf6-4db3-b616-fbd5d50595b1
begin
	signalpower = zeros(1201,nconWTG,15)
	signalfreq  = zeros(1201,nconWTG,15)
	totalmember = zeros(1,nconWTG)
	for icon = 1 : nconWTG
		if (wtgscheme == "DGW") || (wtgscheme == "LGW")
			fnc = "$wtgscheme-$expname-$(dampingstrprnt(configWTG[icon])).nc"
		else
			fnc = "$wtgscheme-$expname-$(relaxscalestrprnt(configWTG[icon])).nc"
		end
		nmem = 0
		if isfile(datadir("precipitation",fnc))
			ds_dgwprcp = NCDataset(datadir("precipitation",fnc))
			prcp  = ds_dgwprcp["precipitation"][:,:] / 24
			close(ds_dgwprcp)
			for imem = 1 : 15
				if sum(.!isnan.(prcp[:,imem])) == 6000
					nmem += 1
					prcpii = prcp[(end-2399):end,imem] .- mean(prcp[(end-2399):end,imem])
					pdg = periodogram(prcpii,fs=24)
					signalpower[:,icon,imem] .= pdg.power
					signalfreq[:,icon,imem]  .= pdg.freq
				end
			end
		end
		totalmember[icon] = nmem
	end
	signalpower = dropdims(sum(signalpower,dims=3),dims=3)
	signalpower = signalpower ./ totalmember
	signalfreq  = dropdims(sum(signalfreq ,dims=3),dims=3)
	signalfreq  = signalfreq  ./ totalmember
	signalfreq  = dropdims(mean(signalfreq,dims=2),dims=2)
	md"Doing power spectrum ..."
end

# ╔═╡ d28f2438-8b19-4763-a31d-4cd2feb30ace
begin
	pplt.close(); f2,a2 = pplt.subplots(aspect=4,axwidth=5)

	if (wtgscheme == "DGW") || (wtgscheme == "LGW")
		wtglabel = L"$a_m$ / day$^{-1}$"
	else
		wtglabel = L"$\tau$ / hr"
	end
	
	c2 = a2[1].pcolormesh(1 ./signalfreq[2:end],configWTG,log10.(signalpower')[:,2:end],levels=-5:0.5:-0,extend="both")
	a2[1].format(
		xscale="log",xlim=(0.1,10),
		yscale="log",ylim=(minimum(configWTG),maximum(configWTG)),
		suptitle="$wtgscheme | $expname",
		xlabel=L"Frequency / day$^{-1}$",ylabel=wtglabel
	)

	f2.colorbar(c2,locator=-5:-0)
	f2.savefig(
		plotsdir("03-timeseries-$wtgscheme-$expname.png"),
		transparent=false,dpi=400
	)
	load(plotsdir("03-timeseries-$wtgscheme-$expname.png"))
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
# ╟─3f880613-b30e-417c-99e0-d78ea2302891
# ╟─0a74e728-1cf6-4db3-b616-fbd5d50595b1
# ╟─d28f2438-8b19-4763-a31d-4cd2feb30ace
