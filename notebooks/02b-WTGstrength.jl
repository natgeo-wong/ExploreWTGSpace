### A Pluto.jl notebook ###
# v0.19.8

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
# 2b. Momentum Damping Strength

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
	configDGW = [
		"damping001","damping002","damping004","damping008",
		"damping011","damping016","damping023","damping032",
		"damping045","damping064","damping090","damping128",
		"damping256","damping512",
	]
	configWTG = [
		"relaxscale02d0","relaxscale03d2","relaxscale04d0","relaxscale05d0",
		"relaxscale05d9","relaxscale07d1","relaxscale08d4","relaxscale10d0",
		"relaxscale11d9","relaxscale14d1","relaxscale16d8","relaxscale20d0",
		"relaxscale25d1","relaxscale31d6","relaxscale50d0",
	]
	nconDGW = length(configDGW)
	nconWTG = length(configWTG)
	bluesDGW = pplt.get_colors("Blues",(nconDGW+2))
	bluesWTG = pplt.get_colors("Blues",(nconWTG+2))
	grnsDGW  = pplt.get_colors("Teal",(nconDGW+2))
	grnsWTG  = pplt.get_colors("Teal",(nconWTG+2))
	brwnsDGW = pplt.get_colors("Brown",(nconDGW+2))
	brwnsWTG = pplt.get_colors("Brown",(nconWTG+2))
	lgd = Dict("frame"=>false,"ncols"=>4)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
# ╠═╡ show_logs = false
begin
	pplt.close()
	fts,ats = pplt.subplots(
		ncols=6,nrows=2,aspect=0.4,axwidth=1.2,
		sharex=0,sharey=0,wspace=1.5
	)
	
	for imem = 1 : 10
		fnc = outstatname(expname,"",true,false,false,false,true,imem)
		if isfile(fnc)
			_,p,t = retrievedims_fnc(fnc); t = t .- 80
			pr = retrievevar_fnc("PREC",fnc)./24
			pa = retrievevar_fnc("AREAPREC",fnc)
			pw = retrievevar_fnc("PW",fnc)
			ta = retrievevar_fnc("TABS",fnc)
			qv = retrievevar_fnc("QV",fnc)
			ol = retrievevar_fnc("LWNT",fnc)
			rh = calcrh(qv,ta,p)
			sw = calcswp(rh,qv,p)
			cr = pw ./ sw / 10
			pra = pr./pa; pra[isnan.(pra)] .= 0; pra[pra.==Inf] .= 0
			ats[1].plot([1,1]*mean(pr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[2].plot([1,1]*mean(pa[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[3].plot([1,1]*mean(pra[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[4].plot([1,1]*mean(pw[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[5].plot([1,1]*mean(cr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[6].plot([1,1]*mean(ol[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[7].plot([1,1]*mean(pr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[8].plot([1,1]*mean(pa[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[9].plot([1,1]*mean(pra[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[10].plot([1,1]*mean(pw[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[11].plot([1,1]*mean(cr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[12].plot([1,1]*mean(ol[(end-499):end]),[0,2000],c="grey",lw=1)
		end
	end
	
	for ic in 1 : nconDGW
		config = configDGW[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 15; imem += 1
			fnc = outstatname(expname,configDGW[ic],false,true,false,false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims_fnc(fnc); t = t .- 80
				pr = retrievevar_fnc("PREC",fnc)./24
				pa = retrievevar_fnc("AREAPREC",fnc)
				pw = retrievevar_fnc("PW",fnc)
				ra = pr./pa; ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
				ta = retrievevar_fnc("TABS",fnc)
				qv = retrievevar_fnc("QV",fnc)
				ol = retrievevar_fnc("LWNT",fnc)
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
					ats[6].plot(mean(ol[(end-99):end]),config,marker=".",c="k",ms=4)
				else
					ats[1].scatter(mean(pr[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[2].scatter(mean(pa[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[3].scatter(mean(ra[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[4].scatter(mean(pw[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[5].scatter(mean(cr[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[6].scatter(mean(ol[(end-99):end]),config,c=clr,alpha=0.2,s=50)
				end
			end
		end
		
	end

	for ic in 1 : nconWTG
		config = configWTG[ic]
		config = replace(config,"relaxscale"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 15; imem += 1
			fnc = outstatname(expname,configWTG[ic],false,false,true,false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims_fnc(fnc); t = t .- 80
				pr = retrievevar_fnc("PREC",fnc)./24
				pa = retrievevar_fnc("AREAPREC",fnc)
				pw = retrievevar_fnc("PW",fnc)
				ra = pr./pa; ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
				ta = retrievevar_fnc("TABS",fnc)
				qv = retrievevar_fnc("QV",fnc)
				ol = retrievevar_fnc("LWNT",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[7].plot(mean(pr[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[8].plot(mean(pa[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[9].plot(mean(ra[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[10].plot(mean(pw[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[11].plot(mean(cr[(end-99):end]),config,marker=".",c="k",ms=4)
					ats[12].plot(mean(ol[(end-99):end]),config,marker=".",c="k",ms=4)
				else
					ats[7].scatter(mean(pr[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[8].scatter(mean(pa[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[9].scatter(mean(ra[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[10].scatter(mean(pw[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[11].scatter(mean(cr[(end-99):end]),config,c=clr,alpha=0.2,s=50)
					ats[12].scatter(mean(ol[(end-99):end]),config,c=clr,alpha=0.2,s=50)
				end
			end
		end
		
	end
	
	ats[1].format(
		ylabel=L"$a_m$ / day$^{-1}$",
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,10),ultitle="(a)",
		ltitle=L"(1) DGW Implementation | Sensitivity to $a_m$ | " * "$(expname)",
	)
	
	ats[2].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,1),ultitle="(b)",
	)
	
	ats[3].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>1),
		xlim=(0,10),ultitle="(c)",
	)
	
	ats[4].format(xlim=(0,75),ultitle="(d)")
	ats[5].format(xlim=(0,100),ultitle="(e)")
	ats[6].format(xlim=(350,100),ultitle="(f)")

	for ii in 1 : 6
		ats[ii].format(ylim=(0.5,2000),yscale="log")
	end
	for ii in 2 : 6
		ats[ii].format(yticklabels=["","","",""],ytickminor=10:10:100)
	end

	ats[7].format(
		ylabel=L"$\tau$ / hr",
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,10),xlabel=L"Domain $P$ / mm hr$^{-1}$",ultitle="(a)",
		ltitle=L"(2) WTG Implementation | Sensitivity to $a_m$ | " * "$(expname)",
	)
	
	ats[8].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
		xlim=(0,1),xlabel="Rain Area Fraction",ultitle="(b)",
	)
	
	ats[9].format(
		xscale="symlog",xscale_kw=Dict("linthresh"=>1),
		xlim=(0,10),xlabel=L"Rain Area $P$ / mm hr$^{-1}$",ultitle="(c)",
	)
	
	ats[10].format(xlim=(0,75),xlabel="PWV / mm",ultitle="(d)")
	ats[11].format(xlim=(0,100),xlabel="CRH / %",ultitle="(e)")
	ats[12].format(xlim=(350,100),xlabel=L"OLR / W m$^{-2}$",ultitle="(f)")

	for ii in 7 : 12
		ats[ii].format(ylim=(2,50),yscale="log")
	end
	for ii in 8 : 12
		ats[ii].format(yticklabels=["","","",""],ytickminor=10:10:100)
	end
	
	for ax in ats
		ax.format(lrtitle="Wet",lltitle="Dry")
	end
	
	fts.savefig(plotsdir(
		"02b-bifurcation-$(expname).png"),
		transparent=true,dpi=400
	)
	load(plotsdir("02b-bifurcation-$(expname).png"))
end

# ╔═╡ 9cf4fa56-91a8-11eb-2710-955eefd10142
md"
We see the following:
* As $a_m$ increases, both the wet and dry model regimes converge back into the initial RCE state.
* When $a_m$ becomes too small, all model states collapse into a dry regime.
* As $a_m$ decreases, the change in rain-area fraction in the domain, rather than rain-rate in the rain-area, which is responsible for changes in the domain-averaged precipitation.

In conclusion, the overall transition from RCE to WTG forcing (decreasing $a_m$) is as follows:
1. Domain becomes more moist / heavier precipitation / wetter
2. Eventually, a dry regime separates out
3. Wet regime begins to show drier characteristics
4. Eventually, wet regime crosses a threshold and model fully enters into a dry-regime.

We note that the simulations that end up in the wet regime are more numerous than those that end up in the dry regime.  Overall, assuming complete randomness this seems to indicate that the a wet regime is favoured.

I have decided on using precipitable water as the prognostic variable for determining if the model is in a dry or wet regime.
"

# ╔═╡ 364a1ce8-91ba-11eb-29a8-b948110e6125
md"
### B. Exploring some Variables
"

# ╔═╡ 489b5bea-91b4-11eb-358b-3fe61c900511
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
		
		while imem < 15; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
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
		"02b-seb-$(expname).png"),
		transparent=false,dpi=300
	)
	load(plotsdir("02b-seb-$(expname).png"))
end

# ╔═╡ e967eb5c-91c4-11eb-3066-05ccaa40bd11
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
		
		while imem < 15; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
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
				# a3D[2].scatter(dropdims(tabi,dims=2),p,c="k")
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
		"02b-vertprofiles-$(expname).png"),
		transparent=false,dpi=200
	)
	load(plotsdir("02b-vertprofiles-$(expname).png"))
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
# ╠═55230f4a-7661-11eb-1c37-8b022b95e08e
# ╟─9cf4fa56-91a8-11eb-2710-955eefd10142
# ╟─364a1ce8-91ba-11eb-29a8-b948110e6125
# ╟─489b5bea-91b4-11eb-358b-3fe61c900511
# ╟─e967eb5c-91c4-11eb-3066-05ccaa40bd11
