### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

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
# 2d. A Possible 2-Day Wave Behaviour

We have previously seen in notebook `02c` that the DGW implementation in idealized conditions gives rise to something approximating a 2-Day Wave Behaviour. In this notebook we investigate the timeseries of precipitation and column water vapour and perform a spectral analysis to confirm this point.
"

# ╔═╡ b6892634-9199-11eb-38d5-8dda8da16ed7
md"
### A. Performing Spectral Analysis

Text
"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configDGW = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500,1000]
	configWTG = [
		0.01,0.02,0.02*sqrt(2.5),0.05,0.05*sqrt(2),
		0.1,0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),
		1,sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2),100
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

# ╔═╡ 56cb6812-df09-411c-9d2a-b49b7dcd6d08
begin
	signalpower_DGW = zeros(401,nconDGW,15)
	signalfreq_DGW  = zeros(401,nconDGW,15)
	totalmember_DGW = zeros(1,nconDGW)
	for icon = 1 : nconDGW
		nmem = 0
		for imem = 1 : 15
			config = "damping$(dampingstrprnt(configDGW[icon]))"
			fnc = outstatname("DGW","S1284km300V64",config,false,true,imem)
			if isfile(fnc)
				prcp = retrievevar_fnc("PREC",fnc)
				if sum(.!isnan.(prcp)) == 2000
					nmem += 1
					prcp = prcp[(end-799):end] .- mean(prcp[(end-799):end])
					pdg = periodogram(prcp,fs=8)
					signalpower_DGW[:,icon,imem] .= pdg.power
					signalfreq_DGW[:,icon,imem]  .= pdg.freq
				end
			end
		end
		totalmember_DGW[icon] = nmem
	end
	signalpower_DGW = dropdims(sum(signalpower_DGW,dims=3),dims=3)
	signalpower_DGW = signalpower_DGW ./ totalmember_DGW
	signalfreq_DGW  = dropdims(sum(signalfreq_DGW ,dims=3),dims=3)
	signalfreq_DGW  = signalfreq_DGW  ./ totalmember_DGW
	signalfreq_DGW  = dropdims(mean(signalfreq_DGW,dims=2),dims=2)
end

# ╔═╡ d50872fc-4e67-4c6b-9d3f-4f478fbfb1f7
begin
	signalpower_TGR = zeros(401,nconWTG,15)
	signalfreq_TGR  = zeros(401,nconWTG,15)
	totalmember_TGR = zeros(1,nconWTG)
	for icon = 1 : nconWTG
		nmem = 0
		for imem = 1 : 15
			config = relaxscalestrprnt(configWTG[icon])
			fnc = outstatname("WTG","S1284km300V64",config,false,true,imem)
			if isfile(fnc)
				prcp = retrievevar_fnc("PREC",fnc)
				if sum(.!isnan.(prcp)) == 2000
					nmem += 1
					prcp = prcp[(end-799):end] .- mean(prcp[(end-799):end])
					pdg = periodogram(prcp,fs=8)
					signalpower_TGR[:,icon,imem] .= pdg.power
					signalfreq_TGR[:,icon,imem]  .= pdg.freq
				end
			end
		end
		totalmember_TGR[icon] = nmem
	end
	signalpower_TGR = dropdims(sum(signalpower_TGR,dims=3),dims=3)
	signalpower_TGR = signalpower_TGR ./ totalmember_TGR
	signalfreq_TGR  = dropdims(sum(signalfreq_TGR ,dims=3),dims=3)
	signalfreq_TGR  = signalfreq_TGR  ./ totalmember_TGR
	signalfreq_TGR  = dropdims(mean(signalfreq_TGR,dims=2),dims=2)
end

# ╔═╡ b8c4ecaf-c1b4-4a6f-b32c-1424cbcbe2f5
begin
	pplt.close(); fig,axs = pplt.subplots(ncols=2,aspect=2,axwidth=3,sharey=0)

	lvls = 10:10:100
	
	c = axs[1].contourf(
		1 ./signalfreq_DGW[2:end],configDGW,
		(signalpower_DGW)'[:,2:end],
		levels=lvls,extend="both",cmap="fire"
	)
	axs[1].format(
		ylabel=L"$a_m$ / day$^{-1}$",ylim=(0.01,100),
		ultitle="(c) Power Spectral Density (DGW)"
	)
	
	c = axs[2].contourf(
		1 ./signalfreq_TGR[2:end],configWTG,
		(signalpower_TGR)'[:,2:end],
		levels=lvls,extend="both",cmap="fire"
	)
	axs[2].format(
		ylabel=L"$\tau$ / hr",ylim=(0.01,100),
		ultitle="(d) Power Spectral Density (TGR)",ytickloc="right"
	)

	for ax in axs
		ax.format(
			xscale="log",xlim=(0.2,50),xlabel="Period / Days",
			yscale="log",yformatter="log"
		)
	end

	fig.colorbar(c,length=0.75,label="Power / dB")
	fig.savefig(plotsdir("02d-2DayWave.png"),transparent=false,dpi=400)
	load(plotsdir("02d-2DayWave.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─56cb6812-df09-411c-9d2a-b49b7dcd6d08
# ╟─d50872fc-4e67-4c6b-9d3f-4f478fbfb1f7
# ╟─b8c4ecaf-c1b4-4a6f-b32c-1424cbcbe2f5
