### A Pluto.jl notebook ###
# v0.19.8

using Markdown
using InteractiveUtils

# ╔═╡ b922b1af-31e6-4ea1-b30b-c5f03f882857
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 181a3e2a-126b-4392-a2de-b2dbd68c9d15
begin
	@quickactivate "ExploreWTGSpace"
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))

	md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ d41d44c4-e9d0-11ec-0b8d-0d5edc09dbd2
md"
# 03b. Creating a WTG Schematic
"

# ╔═╡ 89bb219c-0fab-4126-a529-20c487023a7e
expname = "P1282km300V64"

# ╔═╡ 8c75f5c0-37d8-45e9-b361-709e430fed25
config = "damping016"

# ╔═╡ a29eda09-ac40-4131-8d86-e2ba55b50976
md"
### A. Loading the RCE, Large-Scale and WTG Data
"

# ╔═╡ 5eda6b78-22dc-4554-9a5d-ac213a9fa172
begin
	rh_RCE = zeros(64,2000)
	p_RCE  = zeros(64)
	for imem = 1 : 10
		_,p,t = retrievedims(expname,isensemble=true,member=imem); t = t .- 80
		ta = retrievevar("TABS",expname,isensemble=true,member=imem)
		qv = retrievevar("QV",expname,isensemble=true,member=imem) / 1000
		rh_RCE[:,:] += calcrh(qv,ta,p) * 100
		p_RCE[:]    += p
	end
	rh_RCE = mean(rh_RCE[:,(end-499):end],dims=2) / 10
	p_RCE  = p_RCE / 10
	md"Loading relative humidity profile for RCE simulations ..."
end

# ╔═╡ 0a73df1e-4244-47a0-bafb-b79c0a2d3357
begin
	rh_dry = zeros(64,500)
	w_dry  = zeros(64,500)
	p_dry  = zeros(64)
	for imem = 6 : 10
		_,p,_ = retrievedims(expname,config,isDGW=true,isensemble=true,member=imem)
		ta = retrievevar("TABS",expname,config,isDGW=true,isensemble=true,member=imem)
		qv = retrievevar("QV",expname,config,isDGW=true,isensemble=true,member=imem)
		w  = retrievevar("WWTG",expname,config,isDGW=true,isensemble=true,member=imem)
		rh_dry[:,:] += calcrh(qv,ta,p) / 10
		p_dry[:]    += p
		w_dry[:,:]  += w
	end
	rh_dry = mean(rh_dry[:,(end-99):end],dims=2) / 5
	w_dry  = mean(w_dry[:,(end-99):end],dims=2)  / 5
	p_dry  = p_dry / 5
	md"Loading relative humidity profile for dry WTG simulations ..."
end

# ╔═╡ 91e9827a-db28-4ba8-a2dc-a76e71e21993
begin
	rh_wet = zeros(64,500)
	w_wet  = zeros(64,500)
	p_wet  = zeros(64)
	for imem = 11 : 15
		_,p,_ = retrievedims(expname,config,isDGW=true,isensemble=true,member=imem)
		ta = retrievevar("TABS",expname,config,isDGW=true,isensemble=true,member=imem)
		qv = retrievevar("QV",expname,config,isDGW=true,isensemble=true,member=imem)
		w  = retrievevar("WWTG",expname,config,isDGW=true,isensemble=true,member=imem)
		rh_wet[:,:] += calcrh(qv,ta,p) / 10
		p_wet[:]    += p
		w_wet[:,:]  += w
	end
	rh_wet = mean(rh_wet[:,(end-99):end],dims=2) / 5
	w_wet  = mean(w_wet[:,(end-99):end],dims=2)  / 5
	p_wet  = p_wet / 5
	md"Loading relative humidity profile for wet WTG simulations ..."
end

# ╔═╡ c76f9180-d477-4ec4-84f5-ac2d6adb95b3
begin
	fnc_2D = "RCE_SelfAggregation-Minute5_64_000000"
	olr = zeros(2048,64)
	x   = zeros(2048)
	y   = zeros(64)
	nhr = 12
	for it = (288 - nhr + 1) : 288
		fnc = datadir("LSD","OUT_2D","$(fnc_2D)$(it*10).2Dbin_1.nc")
		d2D = NCDataset(fnc)
		x[:] += d2D["x"][:] / 1000
		y[:] += d2D["y"][:] / 1000
		olr[:,:] += d2D["LWNT"][:]
		close(d2D)
	end
	x = x / nhr
	y = y / nhr
	olr = olr / nhr
	md"Loading Large-Scale data ..."
end

# ╔═╡ 16fb08d4-7bbc-4422-8755-530528e7a8f1
begin
	dSS = NCDataset(datadir("SSD","OUT_2D","SAMTEST_SelfAggregation-NoMPI_32_0000288000.2Dbin_1.nc"))
	oss = dSS["LWNT"][:,:,1]
	xss = dSS["x"][:] / 1000
	yss = dSS["y"][:] / 1000
	close(dSS)
end

# ╔═╡ 41c12fde-938b-40f3-acf6-fecbae34e21c
begin
	arr = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
		[5,7,7,7,7,6,11,11,11,2,4,4,4,4,3,12,12,12,8,10,10,10,10,9]
	]
	pplt.close()
	fig,axs = pplt.subplots(arr,aspect=6,axwidth=6,wspace=0.2,sharex=0,sharey=0)

	c1 = axs[1].pcolormesh(
		x[669:1052],y,olr[669:1052,:]',
		levels=120:10:270,extend="both",cmap="Blues"
	)
	axs[1].plot(x[[973,1004,1004,973,973]],y[[60,60,29,29,60]],lw=3,c="k")
	axs[1].plot(x[[713,744,744,713,713]],y[[38,38,7,7,38]],lw=3,c="w")
	axs[1].text(1432,54,"(c)",c="w",fontweight="bold")
	axs[1].text(1952,100,"(d)",c="k",fontweight="bold")
	axs[1].format(
		grid="false",ylim=(0,128),ylocator=0:32:128,
		xlabel="X / km",ylabel="Y / km",ltitle="(a) Large-Domain RCE"
	)

	axs[2].plot([0,0],[1000,25])
	axs[2].format(
		xlim=(-1,1),ylim=(1000,25),yscale="symlog",xscale="symlog",
		xscale_kw=Dict("linthresh"=>0.01),xlocator=[0],xticklabels=["0"],
		yloc="neither",xloc="bottom",xlabel=L"$w_{wtg}$"
	)

	axs[3].plot(rh_RCE,p_RCE)
	axs[3].format(
		xlim=(-20,120),ylim=(1000,25),yscale="log",
		xloc="bottom",yloc="neither",xlabel="r / %"
	)

	axs[4].contourf(
		xss[65:128],yss,oss[65:128,:]',
		levels=120:10:270,extend="both",cmap="Blues"
	)
	axs[4].format(
		xtickloc="neither",ytickloc="neither",
		grid="false",title="(b) Small-Domain RCE"
	)

	axs[5].plot(w_dry,p_dry)
	axs[5].format(
		xlim=(-1,1),ylim=(1000,25),yscale="symlog",xscale="symlog",
		xscale_kw=Dict("linthresh"=>0.01),xlocator=[0],xticklabels=["0"],
		yloc="left",xloc="bottom",xlabel=L"$w_{wtg}$",ylabel="p / hPa"
	)
	
	axs[6].plot(rh_RCE,p_RCE,c="gray3",linestyle="-.")
	axs[6].plot(rh_dry,p_dry)
	axs[6].format(
		xlim=(-20,120),ylim=(1000,25),yscale="log",
		xloc="bottom",yloc="neither",xlabel="r / %"
	)
	
	axs[7].contourf(
		x[713:744],y[7:38],olr[713:744,7:38]',
		levels=120:10:270,extend="both",cmap="Blues"
	)
	axs[7].format(
		xtickloc="neither",ytickloc="neither",
		grid="false",title="(c) Dry Regime"
	)
	
	axs[8].plot(w_wet,p_wet)
	axs[8].format(
		xlim=(-1,1),ylim=(1000,25),yscale="symlog",xscale="symlog",
		xscale_kw=Dict("linthresh"=>0.01),xlocator=[0],xticklabels=["0"],
		yloc="neither",xloc="bottom",xlabel=L"$w_{wtg}$"
	)

	axs[9].plot(rh_RCE,p_RCE,c="gray3",linestyle="-.")
	axs[9].plot(rh_wet,p_wet)
	axs[9].format(
		xlim=(-20,120),ylim=(1000,25),yscale="log",
		xloc="bottom",yloc="right",xlabel="r / %",ylocator=[]
	)

	axs[10].contourf(
		x[973:1004],y[29:60],olr[973:1004,29:60]',
		levels=120:10:270,extend="both",cmap="Blues"
	)
	axs[10].format(
		xtickloc="neither",ytickloc="neither",
		grid="false",title="(d) Wet Regime"
	)

	for ii in [11,12]
		for addii = 1.55 : 0.15 : 2.95
			axs[ii].plot(-0.8:0.01:1,sin.((-0.8:0.01:1)*pi*5.5)*0.02 .+addii,c="gray4")
		end
		axs[ii].format(xloc="neither",yloc="neither",xticks=[],yticks=[],xlim=(-1,1))
	end

	fig.colorbar(c1,label=L"Outgoing Longwave Radiation / W m$^{-2}$")
	fig.savefig(plotsdir("03b-WTGSchematic.png"),transparent=true,dpi=600)
	load(plotsdir("03b-WTGSchematic.png"))
end

# ╔═╡ Cell order:
# ╟─d41d44c4-e9d0-11ec-0b8d-0d5edc09dbd2
# ╟─b922b1af-31e6-4ea1-b30b-c5f03f882857
# ╟─181a3e2a-126b-4392-a2de-b2dbd68c9d15
# ╠═89bb219c-0fab-4126-a529-20c487023a7e
# ╠═8c75f5c0-37d8-45e9-b361-709e430fed25
# ╟─a29eda09-ac40-4131-8d86-e2ba55b50976
# ╟─5eda6b78-22dc-4554-9a5d-ac213a9fa172
# ╟─0a73df1e-4244-47a0-bafb-b79c0a2d3357
# ╟─91e9827a-db28-4ba8-a2dc-a76e71e21993
# ╟─c76f9180-d477-4ec4-84f5-ac2d6adb95b3
# ╟─16fb08d4-7bbc-4422-8755-530528e7a8f1
# ╟─41c12fde-938b-40f3-acf6-fecbae34e21c
