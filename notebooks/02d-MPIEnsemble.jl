### A Pluto.jl notebook ###
# v0.19.0

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
# 2d. Using the MPI Ensemble in Perturbations

There is a thought that maybe, especially in the simulations with static radiation and fixed surface fluxes, that the perturbations away from the RCE states at certain $a_m$ values is due to random noise that is quite evident in large domains.  As a result, we have introduced MPI Ensembles, where individual subdomains are run on individual processes, with only the large-scale forcing being based on domain-wide averages.  This would help to reduce random noise in our perturbations.
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
	"S" => "Bulk-surface Fluxes (S)",
])

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
begin	
	expname = "$(prefix)1284km300V64"
	
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
	blues = pplt.get_colors("Blues",(ncon+2))
	grns  = pplt.get_colors("Teal",(ncon+2))
	brwns = pplt.get_colors("Brown",(ncon+2))
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
		
		while imem < 3; imem += 1
			fnc = outstatname(expname,configvec[ic],true,true,imem)
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

				if imem == 2
					clr = "yellow7"
				elseif imem == 3
					clr = "blue5"
				end

				if imem == 1
					ats[1].plot(mean(pr[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[2].plot(mean(pa[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[3].plot(mean(ra[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[4].plot(mean(pw[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[5].plot(mean(cr[(end-99):end]),config,marker=".",c="k",ms=4)
				else
					ats[1].scatter(mean(pr[(end-99):end]),config,c=clr,alpha=0.5,s=50)
					ats[2].scatter(mean(pa[(end-99):end]),config,c=clr,alpha=0.5,s=50)
					ats[3].scatter(mean(ra[(end-99):end]),config,c=clr,alpha=0.5,s=50)
					ats[4].scatter(mean(pw[(end-99):end]),config,c=clr,alpha=0.5,s=50)
					ats[5].scatter(mean(cr[(end-99):end]),config,c=clr,alpha=0.5,s=50)
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
		"02d-MPIEnsemble-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("02d-MPIEnsemble-$(expname).png"))
end

# ╔═╡ 13d4b942-0714-494a-bb2c-362ad13abd9e
begin
	pplt.close()
	feb,aeb = pplt.subplots(ncols=5,aspect=0.4,axwidth=1.2,sharex=0); pw = zeros(10)
	
	for imem = 1 : 10
		fnc = outstatname("Control",expname,false,true,imem)
		if isfile(fnc)
			_,_,t = retrievedims(fnc); t = t .- 80
			sw = retrievevar("SWNS",fnc)
			lw = -retrievevar("LWNS",fnc)
			sh = -retrievevar("SHF",fnc)
			lh = -retrievevar("LHF",fnc)
			sb = sw .+ lw .+ sh .+ lh
			pw[imem] = mean(retrievevar("PREC",fnc)[(end-499):end])
			aeb[1].plot([1,1]*mean(sw[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[2].plot([1,1]*mean(lw[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[3].plot([1,1]*mean(sh[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[4].plot([1,1]*mean(lh[(end-499):end]),[0,2000],c="grey",lw=1)
			aeb[5].plot([1,1]*mean(sb[(end-499):end]),[0,2000],c="grey",lw=1)
		end
	end
	
	pw = mean(pw)
	
	for ic in 1 : ncon
		icon = configvec[ic]
		icon = replace(icon,"damping"=>"")
		icon = replace(icon,"d"=>".")
		icon = parse(Float64,icon)
		imem = 0
		
		while imem < 5; imem += 1
			fnc = outstatname(expname,configvec[ic],true,true,imem)
			if isfile(fnc)
				_,_,t = retrievedims(fnc); t = t .- 80
				sw = mean(retrievevar("SWNS",fnc)[(end-99):end])
				lw = -mean(retrievevar("LWNS",fnc)[(end-99):end])
				sh = -mean(retrievevar("SHF",fnc)[(end-99):end])
				lh = -mean(retrievevar("LHF",fnc)[(end-99):end])
				sb = sw .+ lw .+ sh .+ lh
				pwi = mean(retrievevar("PREC",fnc)[(end-99):end])
				if pwi < (0.9 * pw)
					aeb[1].scatter(sw,icon,color="yellow7",alpha=0.5,lw=0)
					aeb[2].scatter(lw,icon,color="yellow7",alpha=0.5,lw=0)
					aeb[3].scatter(sh,icon,color="yellow7",alpha=0.5,lw=0)
					aeb[4].scatter(lh,icon,color="yellow7",alpha=0.5,lw=0)
					aeb[5].scatter(sb,icon,color="yellow7",alpha=0.5,lw=0)
				elseif pwi > (1.1 * pw)
					aeb[1].scatter(sw,icon,color="blue5",alpha=0.5,lw=0)
					aeb[2].scatter(lw,icon,color="blue5",alpha=0.5,lw=0)
					aeb[3].scatter(sh,icon,color="blue5",alpha=0.5,lw=0)
					aeb[4].scatter(lh,icon,color="blue5",alpha=0.5,lw=0)
					aeb[5].scatter(sb,icon,color="blue5",alpha=0.5,lw=0)
					
				else
					aeb[1].scatter(sw,icon,color="teal5",alpha=0.5,lw=0)
					aeb[2].scatter(lw,icon,color="teal5",alpha=0.5,lw=0)
					aeb[3].scatter(sh,icon,color="teal5",alpha=0.5,lw=0)
					aeb[4].scatter(lh,icon,color="teal5",alpha=0.5,lw=0)
					aeb[5].scatter(sb,icon,color="teal5",alpha=0.5,lw=0)
				end
			end
		end
		
	end
	
	aeb[1].format(
		ylim=(0.5,2000),ylabel=L"$a_m$ / day$^{-1}$",yscale="log",
		xlim=(0,400),xlabel="Net Shortwave",
		suptitle=L"Energy Balance / W m$^{-2}$ | " * "$(expname)",
		ultitle="(a)"
	)
	
	aeb[2].format(xlim=(-150,0),xlabel="Net Longwave",ultitle="(b)")
	aeb[3].format(xlim=(-75,0),xlabel="Sensible Heat",ultitle="(c)")
	aeb[4].format(xlim=(-300,0),xlabel="Latent Heat",ultitle="(d)")
	aeb[5].format(xlim=(-350,350),xlabel="Surface Balance",ultitle="(e)")
	
	feb.savefig(plotsdir(
		"02d-mpiseb-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("02d-mpiseb-$(expname).png"))
end

# ╔═╡ 32e9b932-384b-49c5-bdab-62f492e86500
begin
	pplt.close()
	f3D,a3D = pplt.subplots(ncols=5,aspect=0.4,axwidth=1.2,sharex=0)
	clc = zeros(64)
	tab = zeros(64)
	swh = zeros(64)
	
	for ic in 1 : ncon
		icon = configvec[ic]
		icon = replace(icon,"damping"=>"")
		icon = replace(icon,"d"=>".")
		icon = parse(Float64,icon)
		imem = 0
		
		while imem < 3; imem += 1
			fnc = outstatname(expname,configvec[ic],true,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims(fnc); t = t .- 80
				clci = mean(retrievevar("CLD",fnc)[:,(end-99):end],dims=2)*100
				tabi = mean(retrievevar("TABS",fnc)[:,(end-99):end],dims=2)
				tabo = mean(retrievevar("TABSOBS",fnc)[:,(end-99):end],dims=2)
				qvi  = mean(retrievevar("QV",fnc)[:,(end-99):end],dims=2) / 10
				rhi  = calcrh(qvi,tabi,p)
				wwtg = mean(retrievevar("WWTG",fnc)[:,(end-99):end],dims=2)
				pwi  = mean(retrievevar("PREC",fnc)[(end-99):end])
				if pwi < (0.9 * pw)
					a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[3].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[4].plot(dropdims(wwtg,dims=2),p,lw=1,c=brwns[ic+1])
					a3D[5].plot(dropdims(rhi,dims=2),p,lw=1,c=brwns[ic+1])
				elseif pwi > (1.1 * pw)
					a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=blues[ic+1])
					a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=blues[ic+1])
					a3D[3].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=blues[ic+1])
					a3D[4].plot(dropdims(wwtg,dims=2),p,lw=1,c=blues[ic+1])
					a3D[5].plot(dropdims(rhi,dims=2),p,lw=1,c=blues[ic+1])
				else
					a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=grns[ic+1])
					a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=grns[ic+1])
					a3D[3].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=grns[ic+1])
					a3D[4].plot(dropdims(wwtg,dims=2),p,lw=1,c=grns[ic+1])
					a3D[5].plot(dropdims(rhi,dims=2),p,lw=1,c=grns[ic+1])
				end
			end
		end
		
	end
	
	for imem = 1 : 10
		fnc = outstatname("Control",expname,false,true,imem)
		if isfile(fnc)
			_,p,_ = retrievedims(fnc);
			clc[:] = mean(retrievevar("CLD",fnc)[:,(end-499):end],dims=2) * 100
			tab[:] = mean(retrievevar("TABS",fnc)[:,(end-499):end],dims=2)
			# swh[:] = mean(retrievevar("RADQRSW",fnc)[:,(end-499):end],dims=2) .+ mean(retrievevar("RADQRLW",fnc)[:,(end-499):end],dims=2)
			qv = mean(retrievevar("QV",fnc)[:,(end-99):end],dims=2) / 10
			rh = calcrh(qv,tab,p)
			a3D[1].plot(clc,p,color="k")
			a3D[2].plot(tab,p,color="k")
			a3D[3].plot(tab*0,p,color="k")
			a3D[4].plot(swh,p,color="k")
			a3D[5].plot(dropdims(rh,dims=2),p,color="k")
		end
	end
	
	a3D[1].format(
		ylim=(1000,20),ylabel="Pressure / hPa",yscale="log",
		xlim=(0,100),xlabel="Cloud Fraction / %",
		suptitle="3D Vertical Profiles | $(expname)",ultitle="(a)"
	)
	
	a3D[2].format(xlim=(150,325),xlabel="T / K",ultitle="(b)")
	a3D[3].format(xlim=(-30,30),xlabel=L"T - T$_{obs}$ / K",ultitle="(c)")
	a3D[4].format(xlim=(-2,2),xlabel=L"$w_{WTG}$ / km hr$^{-1}$",xscale="symlog",
	xscale_kw=Dict("linthresh"=>0.01),ultitle="(c)")
	a3D[5].format(xlim=(0,110),xlabel="Relative Humidity / %",ultitle="(d)")
	# a3D[5].format(xlim=(-350,350),xlabel="Surface Balance",)
	
	f3D.savefig(plotsdir(
		"02d-mpivertprofiles-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("02d-mpivertprofiles-$(expname).png"))
end

# ╔═╡ f7571f52-30d9-4fd4-9882-febf4f26cf04
md"
### B. Let's Have a Look at the Time-Series Plots

For the `S1284km300V64` experiment, we see that there is still a noticeable, consistent WTG forcing.  Here, we compare the time-series of non-mpi and mpi-ensemble simulations to see if there is any difference.
"

# ╔═╡ f18c4dfb-e9e4-4162-8274-9f238d6536c3
begin
	pplt.close()
	f4,a4 = pplt.subplots(nrows=1,aspect=3,axwidth=6,sharey=0)

	for ic in 8
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0

		while imem < 15; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims(fnc); t = t .- floor(t[1])
				pr = retrievevar("PREC",fnc)
				pw = retrievevar("PW",fnc)
				ta = retrievevar("TABS",fnc)
				qv = retrievevar("QV",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10
				if imem <= 5
					a4[1].plot(t,cr,color=grns[4],lw=0.5)
				elseif (imem>=6) && (imem<=10)
					a4[1].plot(t,cr,color=brwns[4],lw=0.5)
				elseif (imem>=11)
					a4[1].plot(t,cr,color=blues[4],lw=0.5)
				end
			end
		end

	end

	for ic in 8
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0

		while imem < 15; imem += 1
			fnc = outstatname(expname,configvec[ic],true,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims(fnc); t = t .- floor(t[1])
				pr = retrievevar("PREC",fnc)
				pw = retrievevar("PW",fnc)
				ta = retrievevar("TABS",fnc)
				qv = retrievevar("QV",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10
				if imem == 1
					a4[1].plot(t,cr,color=grns[10],lw=3)
				elseif imem == 2
					a4[1].plot(t,cr,color=brwns[10],lw=3)
				elseif imem == 3
					a4[1].plot(t,cr,color=blues[10],lw=3)
				end
			end
		end

	end

	a4[1].format(
		ylabel="Column Relative Humidity / %",ylim=(0,100),
		xlim=(0,500),suptitle=expname,ultitle="(a)"
	)

	f4.savefig(plotsdir("02d-rce2wtg-$(expname).png"),transparent=false,dpi=300)
	load(plotsdir("02d-rce2wtg-$(expname).png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╟─e5de2fc0-6f10-4ff9-817f-95fa20821b06
# ╟─d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
# ╟─13d4b942-0714-494a-bb2c-362ad13abd9e
# ╟─32e9b932-384b-49c5-bdab-62f492e86500
# ╟─f7571f52-30d9-4fd4-9882-febf4f26cf04
# ╟─f18c4dfb-e9e4-4162-8274-9f238d6536c3
