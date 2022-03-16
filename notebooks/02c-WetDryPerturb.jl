### A Pluto.jl notebook ###
# v0.18.0

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
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 2c. Wet and Dry Perturbations

In this notebook, we show that even very small perturbations in the large-scale reference profile can lead to a bias in the frequency of the wet and dry states obtained in the bifurcation.  However, we also show that despite changes in the frequency at which the states occur, the climatology of these states remains the same.
"

# ╔═╡ b6892634-9199-11eb-38d5-8dda8da16ed7
md"
### A. Breakdown of the Ensemble Members

In our initial simulations, we ran a 5-member ensemble against the large-scale reference profiles found in notebook `01a`.  Now however, we run 2 more 5-member ensembles for each configuration:
1. A \"hot\" perturbation, +0.05 K to temperature of the large-scale profile
2. A \"cold\" perturbation, -0.05 K to temperature of the large-scale profile

We expect the results of the model from the \"hot\" perturbation to favour the dry state, and vice-versa for the \"cold\" perturbation.
"

# ╔═╡ e5de2fc0-6f10-4ff9-817f-95fa20821b06
@bind prefix Select([
	"P" => "Perpetual Insolation (P)",
	"D" => "Diurnal Insolation (D)",
	"T" => "Non-interactive Radiation (T)",
	"S" => "Bulk-surface Fluxes (S)",
])

# ╔═╡ a99febc5-75f1-416a-8d17-2f6ba4ef9fb0
md"Toggle Domain Size $(@bind islarge PlutoUI.Slider(64:64:128,default=128,show_value=true)) km"

# ╔═╡ ae58e5de-2eb9-4e96-b6ed-2f25c0e682b2
md"Toggle Horizontal Resolution: $(@bind hres PlutoUI.Slider(-1:1,default=0))"

# ╔═╡ ab78df44-4f57-447a-80d3-0531f912a9ed
md"Sea Surface Temperature: $(@bind sst PlutoUI.Slider(295:5:305,default=300, show_value=true))"

# ╔═╡ b8a3a33f-34ca-46c9-867a-f88106ef83cf
md"Coarse Vertical Grid? $(@bind iscvg PlutoUI.Slider(0:1))"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
begin
	domsize = @sprintf("%03d",islarge)
	
	if islarge == 64
		res = 1
	else; res = Int(2. ^hres*2)
	end
	
	if iszero(iscvg)
		  vgrd = 64
	else; vgrd = 28
	end
	
	expname = "$(prefix)$(domsize)$(res)km$(sst)V$(vgrd)"
	
md"**Experiment Set:** $expname"
end

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configvec = [
		"damping001","damping002","damping004","damping008",
		"damping011","damping016","damping023","damping032",
		"damping045","damping064","damping090","damping128",
		"damping256","damping512",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	grns  = pplt.Colors("Teal",(ncon+2))
	brwns = pplt.Colors("Brown",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(ncols=5,aspect=0.4,axwidth=1.2,sharex=0)
	
	for imem = 1 : 10
		fnc = outstatname("Control",expname,false,true,imem)
		if isfile(fnc)
			_,p,t = retrievedims(fnc); t = t .- 80
			pr = retrievevar("PREC",fnc)./24
			pa = retrievevar("AREAPREC",fnc)
			pw = retrievevar("PW",fnc)
			ta = retrievevar("TABS",fnc)
			qv = retrievevar("QV",fnc)
			rh = calcrh(qv,ta,p)
			sw = calcswp(rh,qv,p)
			cr = pw ./ sw / 10
			pra = pr./pa; pra[isnan.(pra)] .= 0; pra[pra.==Inf] .= 0
			ats[1].plot([1,1]*mean(pr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[2].plot([1,1]*mean(pa[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[3].plot([1,1]*mean(pra[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[4].plot([1,1]*mean(pw[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[5].plot([1,1]*mean(cr[(end-499):end]),[0,2000],c="grey",lw=1)
		end
	end
	
	for ic in 1 : ncon
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 15; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims(fnc); t = t .- 80
				pr = retrievevar("PREC",fnc)./24
				pa = retrievevar("AREAPREC",fnc)
				pw = retrievevar("PW",fnc)
				ra = pr./pa; ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
				ta = retrievevar("TABS",fnc)
				qv = retrievevar("QV",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[1].plot(mean(pr[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[2].plot(mean(pa[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[3].plot(mean(ra[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[4].plot(mean(pw[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[5].plot(mean(cr[(end-99):end]),config,marker=".",c="k",ms=4)
				else
					ats[1].scatter(mean(pr[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[2].scatter(mean(pa[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[3].scatter(mean(ra[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[4].scatter(mean(pw[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[5].scatter(mean(cr[(end-99):end]),config,c=clr,alpha=0.2,s=50)
				end
			end
		end
		
	end
	
	ats[1].format(
		ylim=(0.5,2000),ylabel=L"$a_m$ / day$^{-1}$",yscale="log",
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,10),xlabel=L"Domain $P$ / mm hr$^{-1}$",
		suptitle=L"Sensitivity to $a_m$ | " * "$(expname)",
		ultitle="(a)"
	)
	
	ats[2].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,1),xlabel="Rain Area Fraction",
		ultitle="(b)"
	)
	
	ats[3].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>1),
		xlim=(0,10),xlabel=L"Rain Area $P$ / mm hr$^{-1}$",
		ultitle="(c)"
	)
	
	ats[4].format(
		xlim=(0,75),xlabel="PWV / mm",
		ultitle="(d)"
	)
	
	ats[5].format(
		xlim=(0,100),xlabel="CRH / %",
		ultitle="(e)"
	)
	
	for ax in ats
		ax.format(lrtitle="Wet",lltitle="Dry")
	end
	
	fts.savefig(plotsdir(
		"02c-wetdryperturb-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("02c-wetdryperturb-$(expname).png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╟─e5de2fc0-6f10-4ff9-817f-95fa20821b06
# ╟─a99febc5-75f1-416a-8d17-2f6ba4ef9fb0
# ╟─ae58e5de-2eb9-4e96-b6ed-2f25c0e682b2
# ╟─ab78df44-4f57-447a-80d3-0531f912a9ed
# ╟─b8a3a33f-34ca-46c9-867a-f88106ef83cf
# ╟─d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
