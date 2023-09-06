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

# ╔═╡ 3be6a428-9555-4a92-9611-eaf8c0f9d3df
@bind schname Select([
	"VDD" => "(VDD) Vertical-Mode Decomposition for DGW",
	"VDT" => "(VDT) Vertical-Mode Decomposition for TGR",
])

# ╔═╡ e5de2fc0-6f10-4ff9-817f-95fa20821b06
@bind prefix Select([
	"S" => "Bulk-surface Fluxes (S)",
])

# ╔═╡ a99febc5-75f1-416a-8d17-2f6ba4ef9fb0
md"Toggle 1st-Mode Scaling $(@bind mode1 PlutoUI.Slider(0:1,default=1,show_value=true))"

# ╔═╡ ae58e5de-2eb9-4e96-b6ed-2f25c0e682b2
md"Toggle 2nd-Mode Scaling $(@bind mode2 PlutoUI.Slider(0:0.1:1,default=0,show_value=true))"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
begin
	if isone(mode1) && (!isone(mode2) && !iszero(mode2)); mode2_fin = 0.25; else; mode2_fin = mode2 end
	if !isone(mode1); mode2_fin = 1 end
	expname = "H$(@sprintf("%4.2f",mode1))F$(@sprintf("%4.2f",mode2_fin))"
    expname = replace(expname,"."=>"d")
	
md"**WTG Scheme:** $schname | **Experiment Set:** $expname"
end

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configDGW = [1,2,5,10,20,50,100,200,500,1000]
	configTGR = [
		1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100
	]
	if schname == "VDD"
		configWTG = configDGW
	else
		configWTG = configTGR
	end
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
		ncols=5,nrows=2,aspect=0.4,axwidth=1,
		sharex=0,sharey=0,wspace=1.5
	)

	pwr = 0
	for imem = 1 : 10
		fnc = outstatname("RCE","S1284km300V64","",false,true,imem)
		if isfile(fnc)
			_,p,t = retrievedims_fnc(fnc); t = t .- 80
			pr = retrievevar_fnc("PREC",fnc)./24
			pa = retrievevar_fnc("AREAPREC",fnc)
			pw = retrievevar_fnc("PW",fnc)
			global pwr += mean(pr[(end-499):end])
			ta = retrievevar_fnc("TABS",fnc)
			qv = retrievevar_fnc("QV",fnc)
			rh = calcrh(qv,ta,p)
			sw = calcswp(rh,qv,p)
			cr = pw ./ sw / 10
			pra = pr./pa; pra[isnan.(pra)] .= 0; pra[pra.==Inf] .= 0
			ats[1].plot([1,1]*mean(pr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[2].plot([1,1]*mean(pa[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[3].plot([1,1]*mean(pra[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[4].plot([1,1]*mean(pw[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[5].plot([1,1]*mean(cr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[6].plot([1,1]*mean(pr[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[7].plot([1,1]*mean(pa[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[8].plot([1,1]*mean(pra[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[9].plot([1,1]*mean(pw[(end-499):end]),[0,2000],c="grey",lw=1)
			ats[10].plot([1,1]*mean(cr[(end-499):end]),[0,2000],c="grey",lw=1)
		end
	end
	pwr = pwr / 10 * 24
	
	for ic in 1 : nconDGW
		config = "damping$(dampingstrprnt(configDGW[ic]))"
		imem = 0
		
		while imem < 15; imem += 1
			fnc = outstatname("DGW","S1284km300V64",config,false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims_fnc(fnc); t = t .- 80
				pr = retrievevar_fnc("PREC",fnc)./24
				pa = retrievevar_fnc("AREAPREC",fnc)
				pw = retrievevar_fnc("PW",fnc)
				ra = pr./pa; ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
				ta = retrievevar_fnc("TABS",fnc)
				qv = retrievevar_fnc("QV",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10

				pr = mean(pr[(end-99):end])
				pa = mean(pa[(end-99):end])
				ra = mean(ra[(end-99):end])
				pw = mean(pw[(end-99):end])
				cr = mean(cr[(end-99):end])

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[1].plot(pr,configDGW[ic],marker=".",c="k",ms=4)
					ats[2].plot(pa,configDGW[ic],marker=".",c="k",ms=4)
					ats[3].plot(ra,configDGW[ic],marker=".",c="k",ms=4)
					ats[4].plot(pw,configDGW[ic],marker=".",c="k",ms=4)
					ats[5].plot(cr,configDGW[ic],marker=".",c="k",ms=4)
				else
					ats[1].scatter(pr,configDGW[ic],c=clr,alpha=0.1,s=50)
					ats[2].scatter(pa,configDGW[ic],c=clr,alpha=0.1,s=50)
					ats[3].scatter(ra,configDGW[ic],c=clr,alpha=0.1,s=50)
					ats[4].scatter(pw,configDGW[ic],c=clr,alpha=0.1,s=50)
					ats[5].scatter(cr,configDGW[ic],c=clr,alpha=0.1,s=50)
				end
			end
		end
		
	end

	for ic in 1 : nconWTG
		if schname == "VDD"
			config = "damping$(dampingstrprnt(configWTG[ic]))"
		else
			config = relaxscalestrprnt(configWTG[ic])
		end
		imem = 0
		
		while imem < 15; imem += 1
			fnc = outstatname(schname,expname,config,false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				_,p,t = retrievedims_fnc(fnc); t = t .- 80
				pr = retrievevar_fnc("PREC",fnc)./24
				pa = retrievevar_fnc("AREAPREC",fnc)
				pw = retrievevar_fnc("PW",fnc)
				ra = pr./pa; ra[isnan.(ra)] .= 0; ra[ra.==Inf] .= 0
				ta = retrievevar_fnc("TABS",fnc)
				qv = retrievevar_fnc("QV",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10

				pr = mean(pr[(end-99):end])
				pa = mean(pa[(end-99):end])
				ra = mean(ra[(end-99):end])
				pw = mean(pw[(end-99):end])
				cr = mean(cr[(end-99):end])
				

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[6].plot(pr,configWTG[ic],marker=".",c="k",ms=4)
					ats[7].plot(pa,configWTG[ic],marker=".",c="k",ms=4)
					ats[8].plot(ra,configWTG[ic],marker=".",c="k",ms=4)
					ats[9].plot(pw,configWTG[ic],marker=".",c="k",ms=4)
					ats[10].plot(cr,configWTG[ic],marker=".",c="k",ms=4)
				else
					ats[6].scatter(pr,configWTG[ic],c=clr,alpha=0.1,s=50)
					ats[7].scatter(pa,configWTG[ic],c=clr,alpha=0.1,s=50)
					ats[8].scatter(ra,configWTG[ic],c=clr,alpha=0.1,s=50)
					ats[9].scatter(pw,configWTG[ic],c=clr,alpha=0.1,s=50)
					ats[10].scatter(cr,configWTG[ic],c=clr,alpha=0.1,s=50)
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

	for ii in 1 : 5
		ats[ii].format(ylim=(0.5,2000),yscale="log")
	end
	for ii in 2 : 5
		ats[ii].format(yticklabels=["","","",""],ytickminor=10:10:100)
	end

	if schname == "VDD"
		ats[6].format(
			ylabel=L"$a_m$ / day$^{-1}$",
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,10),ultitle="(a)",
			ltitle="(2) " * "$(schname)" * L" Implementation | Sensitivity to $a_m$ | " * "$(expname)",
		)
		
		ats[7].format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,1),ultitle="(b)",
		)
		
		ats[8].format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>1),
			xlim=(0,10),ultitle="(c)",
		)
		
		ats[9].format(xlim=(0,75),ultitle="(d)")
		ats[10].format(xlim=(0,100),ultitle="(e)")

		for ii in 6 : 10
			ats[ii].format(ylim=(0.5,2000),yscale="log")
		end
	else
		ats[6].format(
			ylabel=L"$\tau$ / hr",
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,10),xlabel=L"Domain $P$ / mm hr$^{-1}$",ultitle="(a)",
			ltitle="(2) $(schname)" *L" Implementation | Sensitivity to $\tau$ | " * "$(expname)",
		)
		
		ats[7].format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,1),xlabel="Rain Area Fraction",ultitle="(b)",
		)
		
		ats[8].format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>1),
			xlim=(0,10),xlabel=L"Rain Area $P$ / mm hr$^{-1}$",ultitle="(c)",
		)
		
		ats[9].format(xlim=(0,75),xlabel="PWV / mm",ultitle="(d)")
		ats[10].format(xlim=(0,100),xlabel="CRH / %",ultitle="(e)")

		for ii in 6 : 10
			ats[ii].format(ylim=(0.5,200),yscale="log")
		end
	end
	
	for ii in 7 : 10
		ats[ii].format(yticklabels=["","","",""],ytickminor=10:10:100)
	end
	
	for ax in ats
		ax.format(lrtitle="Wet",lltitle="Dry")
	end
	
	fts.savefig(plotsdir(
		"04b-bifurcation-$(expname)-$(schname).png"),
		transparent=false,dpi=400
	)
	load(plotsdir("04b-bifurcation-$(expname)-$(schname).png"))
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

# ╔═╡ 7f71bb3d-fae1-4d44-a612-0d5bfb1dd4e4
md"Create Image? $(@bind createimage PlutoUI.Slider(0:1))"

# ╔═╡ e967eb5c-91c4-11eb-3066-05ccaa40bd11
begin
	if isone(createimage)
		pplt.close()
		f3D,a3D = pplt.subplots(
			ncols=6,nrows=2,aspect=0.4,axwidth=1,
			sharex=0,wspace=1.5,hspace=4
		)
		clc = zeros(64)
		tab = zeros(64)
		swh = zeros(64)
	
		for ic in 1 : nconDGW
			icon = "damping$(dampingstrprnt(configDGW[ic]))"
			imem = 0
			
			while imem < 15; imem += 1
				fnc = outstatname("DGW","S1284km300V64",icon,false,true,imem)
				if isfile(fnc)
					_,p,t = retrievedims_fnc(fnc); t = t .- 80
					clci = mean(retrievevar_fnc("CLD",fnc)[:,(end-99):end],dims=2)*100
					tabi = mean(retrievevar_fnc("TABS",fnc)[:,(end-99):end],dims=2)
					tabo = mean(retrievevar_fnc("TABSOBS",fnc)[:,(end-99):end],dims=2)
					qvi  = mean(retrievevar_fnc("QV",fnc)[:,(end-99):end],dims=2) / 10
					rhi  = calcrh(qvi,tabi,p)
					wwtg = mean(retrievevar_fnc("WWTG",fnc)[:,(end-99):end],dims=2)
					pwi  = mean(retrievevar_fnc("PREC",fnc)[(end-99):end])
					qvtn = mean(retrievevar_fnc("QVTEND",fnc)[:,(end-99):end],dims=2)
					if pwi < (0.9 * pwr)
						a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=brwnsDGW[ic+1])
						a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=brwnsDGW[ic+1])
						a3D[3].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=brwnsDGW[ic+1])
						a3D[4].plot(dropdims(wwtg,dims=2),p,lw=1,c=brwnsDGW[ic+1])
						a3D[5].plot(dropdims(rhi,dims=2),p,lw=1,c=brwnsDGW[ic+1])
						a3D[6].plot(dropdims(qvtn,dims=2),p,lw=1,c=brwnsDGW[ic+1])
					elseif pwi > (1.1 * pwr)
						a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=bluesDGW[ic+1])
						a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=bluesDGW[ic+1])
						a3D[3].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=bluesDGW[ic+1])
						a3D[4].plot(dropdims(wwtg,dims=2),p,lw=1,c=bluesDGW[ic+1])
						a3D[5].plot(dropdims(rhi,dims=2),p,lw=1,c=bluesDGW[ic+1])
						a3D[6].plot(dropdims(qvtn,dims=2),p,lw=1,c=bluesDGW[ic+1])
					else
						a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=grnsDGW[ic+1])
						a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=grnsDGW[ic+1])
						a3D[3].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=grnsDGW[ic+1])
						a3D[4].plot(dropdims(wwtg,dims=2),p,lw=1,c=grnsDGW[ic+1])
						a3D[5].plot(dropdims(rhi,dims=2),p,lw=1,c=grnsDGW[ic+1])
						a3D[6].plot(dropdims(qvtn,dims=2),p,lw=1,c=grnsDGW[ic+1])
					end
					# a3D[2].scatter(dropdims(tabi,dims=2),p,c="k")
				end
			end
			
		end
	
		for ic in 1 : nconWTG
			if schname == "VDD"
				icon = "damping$(dampingstrprnt(configWTG[ic]))"
			else
				icon = relaxscalestrprnt(configWTG[ic])
			end
			imem = 0
			
			while imem < 15; imem += 1
				fnc = outstatname(schname,expname,icon,false,true,imem)
				fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
				if isfile(fnc)
					_,p,t = retrievedims_fnc(fnc); t = t .- 80
					clci = mean(retrievevar_fnc("CLD",fnc)[:,(end-99):end],dims=2)*100
					tabi = mean(retrievevar_fnc("TABS",fnc)[:,(end-99):end],dims=2)
					tabo = mean(retrievevar_fnc("TABSOBS",fnc)[:,(end-99):end],dims=2)
					qvi  = mean(retrievevar_fnc("QV",fnc)[:,(end-99):end],dims=2) / 10
					rhi  = calcrh(qvi,tabi,p)
					if schname == "VDD"
						wwtg = mean(retrievevar_fnc("OWTG",fnc)[:,(end-99):end],dims=2)
					else
						wwtg = mean(retrievevar_fnc("WWTG",fnc)[:,(end-99):end],dims=2)
					end
					pwi  = mean(retrievevar_fnc("PREC",fnc)[(end-99):end])
					qvtn = mean(retrievevar_fnc("QVTEND",fnc)[:,(end-99):end],dims=2)
					if pwi < (0.9 * pwr)
						a3D[7].plot(dropdims(clci,dims=2),p,lw=1,c=brwnsWTG[ic+1])
						a3D[8].plot(dropdims(tabi,dims=2),p,lw=1,c=brwnsWTG[ic+1])
						a3D[9].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=brwnsWTG[ic+1])
						a3D[10].plot(dropdims(wwtg,dims=2),p,lw=1,c=brwnsWTG[ic+1])
						a3D[11].plot(dropdims(rhi,dims=2),p,lw=1,c=brwnsWTG[ic+1])
						a3D[12].plot(dropdims(qvtn,dims=2),p,lw=1,c=brwnsWTG[ic+1])
					elseif pwi > (1.1 * pwr)
						a3D[7].plot(dropdims(clci,dims=2),p,lw=1,c=bluesWTG[ic+1])
						a3D[8].plot(dropdims(tabi,dims=2),p,lw=1,c=bluesWTG[ic+1])
						a3D[9].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=bluesWTG[ic+1])
						a3D[10].plot(dropdims(wwtg,dims=2),p,lw=1,c=bluesWTG[ic+1])
						a3D[11].plot(dropdims(rhi,dims=2),p,lw=1,c=bluesWTG[ic+1])
						a3D[12].plot(dropdims(qvtn,dims=2),p,lw=1,c=bluesWTG[ic+1])
					else
						a3D[7].plot(dropdims(clci,dims=2),p,lw=1,c=grnsWTG[ic+1])
						a3D[8].plot(dropdims(tabi,dims=2),p,lw=1,c=grnsWTG[ic+1])
						a3D[9].plot(dropdims(tabi.-tabo,dims=2),p,lw=1,c=grnsWTG[ic+1])
						a3D[10].plot(dropdims(wwtg,dims=2),p,lw=1,c=grnsWTG[ic+1])
						a3D[11].plot(dropdims(rhi,dims=2),p,lw=1,c=grnsWTG[ic+1])
						a3D[12].plot(dropdims(qvtn,dims=2),p,lw=1,c=grnsWTG[ic+1])
					end
					# a3D[2].scatter(dropdims(tabi,dims=2),p,c="k")
				end
			end
			
		end
		
		for imem = 1 : 10
			fnc = outstatname("RCE","S1284km300V64","",false,true,imem)
			if isfile(fnc)
				_,p,_ = retrievedims_fnc(fnc);
				clc[:] = mean(retrievevar_fnc("CLD",fnc)[:,(end-499):end],dims=2)*100
				tab[:] = mean(retrievevar_fnc("TABS",fnc)[:,(end-499):end],dims=2)
				qv = mean(retrievevar_fnc("QV",fnc)[:,(end-99):end],dims=2) / 10
				rh = calcrh(qv,tab,p)
				a3D[1].plot(clc,p,color="k")
				a3D[2].plot(tab,p,color="k")
				a3D[3].plot(tab*0,p,color="k")
				a3D[4].plot(swh,p,color="k")
				a3D[5].plot(dropdims(rh,dims=2),p,color="k")
				a3D[6].plot(qv*0,p,color="k")
				a3D[7].plot(clc,p,color="k")
				a3D[8].plot(tab,p,color="k")
				a3D[9].plot(tab*0,p,color="k")
				a3D[10].plot(swh,p,color="k")
				a3D[11].plot(dropdims(rh,dims=2),p,color="k")
				a3D[12].plot(qv*0,p,color="k")
			end
		end
		
		a3D[1].format(
			ylim=(1000,20),ylabel="Pressure / hPa",yscale="log",
			xlim=(0,100),
			ltitle="(1) DGW Implementation | 3D Vertical Profiles | $(expname)",
		)
		
		a3D[2].format(xlim=(150,325),ultitle="(b)")
		a3D[3].format(xlim=(-30,30),ultitle="(c)")
		a3D[4].format(
			xlim=(2,-2),xlabel=L"$w_{WTG}$ / m s$^{-1}$",xscale="symlog",
			xscale_kw=Dict("linthresh"=>0.01),ultitle="(d)"
		)
		a3D[4].format(xlim=(-2,2),xscale="symlog",
		xscale_kw=Dict("linthresh"=>0.01),ultitle="(d)")
		a3D[5].format(xlim=(0,110),ultitle="(e)")
		a3D[6].format(xlim=(-50,50),xscale="symlog",
		xscale_kw=Dict("linthresh"=>1),ultitle="(f)")
	
		a3D[7].format(
			ylim=(1000,20),ylabel="Pressure / hPa",yscale="log",
			xlim=(0,100),xlabel="Cloud Fraction / %",
			ltitle="(2) $(schname) Implementation | 3D Vertical Profiles | $(expname)",
		)
		
		a3D[8].format(xlim=(150,325),xlabel="T / K",ultitle="(b)")
		a3D[9].format(xlim=(-30,30),xlabel=L"T - T$_{obs}$ / K",ultitle="(c)")
		if schname == "VDD"
			a3D[10].format(
				xlim=(2,-2),xlabel=L"$\omega_{WTG}$ / Pa s$^{-1}$",xscale="symlog",
				xscale_kw=Dict("linthresh"=>0.01),ultitle="(d)"
			)
		else
			a3D[10].format(
				xlim=(-2,2),xlabel=L"$w_{WTG}$ / m s$^{-1}$",xscale="symlog",
				xscale_kw=Dict("linthresh"=>0.01),ultitle="(d)"
			)
		end
		a3D[11].format(xlim=(0,110),xlabel="Relative Humidity / %",ultitle="(e)")
		a3D[12].format(xlim=(-50,50),xscale="symlog",
		xscale_kw=Dict("linthresh"=>1),xlabel=L"$\frac{d(H_2O)}{dt}$ / g kg$^{-1}$ day$^{-1}$",ultitle="(f)")
		
		f3D.savefig(plotsdir(
			"04b-vertprofiles-$(expname)-$(schname).png"),
			transparent=false,dpi=200
		)
	end
	load(plotsdir("04b-vertprofiles-$(expname)-$(schname).png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╟─3be6a428-9555-4a92-9611-eaf8c0f9d3df
# ╟─e5de2fc0-6f10-4ff9-817f-95fa20821b06
# ╟─a99febc5-75f1-416a-8d17-2f6ba4ef9fb0
# ╟─ae58e5de-2eb9-4e96-b6ed-2f25c0e682b2
# ╟─d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
# ╟─9cf4fa56-91a8-11eb-2710-955eefd10142
# ╟─364a1ce8-91ba-11eb-29a8-b948110e6125
# ╟─7f71bb3d-fae1-4d44-a612-0d5bfb1dd4e4
# ╟─e967eb5c-91c4-11eb-3066-05ccaa40bd11
