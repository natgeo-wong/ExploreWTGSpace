### A Pluto.jl notebook ###
# v0.19.25

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
	@quickactivate "2023GL104350"
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase

	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")

	include(srcdir("sam.jl"))

	md"Loading modules for the 2023GL104350 paper submission ..."
end

# ╔═╡ d41d44c4-e9d0-11ec-0b8d-0d5edc09dbd2
md"
# Figure 1. A WTG Schematic
"

# ╔═╡ a29eda09-ac40-4131-8d86-e2ba55b50976
md"
### A. Loading the RCE, Large-Scale and WTG Data
"

# ╔═╡ a64e260a-9dd0-491e-b0fc-c4b6fdfa8b6f
function calcT2e(T::Real)

    tb = T - 273.15
    if tb <= 0
        return exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
        return exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end

end

# ╔═╡ 4e4bf212-220d-419e-9e90-fbd524c17881
function calcT2q(T::Real,p::Real)
    e = calcT2e(T)
    return e * 0.621981 / (p - 0.378019 * e)
end

# ╔═╡ 5eda6b78-22dc-4554-9a5d-ac213a9fa172
begin
	ds_rce3D = NCDataset(datadir("3D","RCE-P1282km300V64.nc"))
	p_rce  = ds_rce3D["p"][:] * 100
	qv_rce = ds_rce3D["q"][:,end-999:end,:]
	qv_rce = dropdims(mean(qv_rce,dims=2),dims=2)
	ta_rce = ds_rce3D["t"][:,end-999:end,:]
	ta_rce = dropdims(mean(ta_rce,dims=2),dims=2)
	qs_rce = zeros(size(ta_rce))
	for ien = 1 : 10, ip = 1 : 64
		qs_rce[ip,ien] = calcT2q(ta_rce[ip,ien],p_rce[ip])
	end
	rh_rce = qv_rce ./ qs_rce * 100
	rh_rce = dropdims(mean(rh_rce,dims=2),dims=2)
	close(ds_rce3D)
	md"Loading relative humidity profile for RCE simulations ..."
end

# ╔═╡ 0a73df1e-4244-47a0-bafb-b79c0a2d3357
begin
	ds_dgw3D = NCDataset(datadir("3D","DGW-P1282km300V64-damping010.nc"))
	p_dgw  = ds_dgw3D["p"][:] * 100
	qv_dgw = ds_dgw3D["q"][:,end-999:end,:]
	qv_dgw = dropdims(mean(qv_dgw,dims=2),dims=2)
	ta_dgw = ds_dgw3D["t"][:,end-999:end,:]
	ta_dgw = dropdims(mean(ta_dgw,dims=2),dims=2)
	qs_dgw = zeros(size(ta_dgw))
	for ien = 1 : 15, ip = 1 : 64
		qs_dgw[ip,ien] = calcT2q(ta_dgw[ip,ien],p_dgw[ip])
	end
	rh_dgw = qv_dgw ./ qs_dgw
	rh_dry = dropdims(mean(rh_dgw[:,06:10],dims=2),dims=2) * 100
	rh_wet = dropdims(mean(rh_dgw[:,11:15],dims=2),dims=2) * 100
	close(ds_dgw3D)
	md"Loading relative humidity profile for DGW simulations ..."
end

# ╔═╡ ae968145-c8b4-4e7c-b348-7d20085fac96
begin
	ds_dgww = NCDataset(datadir("wwtg","DGW-P1282km300V64-damping010.nc"))
	w_dgw = ds_dgww["wwtg"][:,end-999:end,:]
	w_dgw = dropdims(mean(w_dgw,dims=2),dims=2)
	w_dry = dropdims(mean(w_dgw[:,06:10],dims=2),dims=2)
	w_wet = dropdims(mean(w_dgw[:,11:15],dims=2),dims=2)
	close(ds_dgww)
	md"Loading vertical velocity profile for DGW simulations ..."
end

# ╔═╡ c76f9180-d477-4ec4-84f5-ac2d6adb95b3
begin
	dLS = NCDataset(datadir("largescale.nc"))
	olr = dLS["olr"][:,:,1]
	x = dLS["x"][:]
	y = dLS["y"][:]
	close(dLS)
	md"Loading Large-Scale data ..."
end

# ╔═╡ 16fb08d4-7bbc-4422-8755-530528e7a8f1
begin
	dSS = NCDataset(datadir("localscale.nc"))
	oss = dSS["olr"][:,:,1]
	xss = dSS["x"][:] / 1000
	yss = dSS["y"][:] / 1000
	close(dSS)
	md"Loading Local-Scale data ..."
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

	axs[3].plot(rh_rce,p_rce/100)
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

	axs[5].plot(w_dry,p_dgw/100)
	axs[5].format(
		xlim=(-1,1),ylim=(1000,25),yscale="symlog",xscale="symlog",
		xscale_kw=Dict("linthresh"=>0.01),xlocator=[0],xticklabels=["0"],
		yloc="left",xloc="bottom",xlabel=L"$w_{wtg}$",ylabel="p / hPa"
	)
	
	axs[6].plot(rh_rce,p_rce/100,c="gray3",linestyle="-.")
	axs[6].plot(rh_dry,p_dgw/100)
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
	
	axs[8].plot(w_wet,p_dgw/100)
	axs[8].format(
		xlim=(-1,1),ylim=(1000,25),yscale="symlog",xscale="symlog",
		xscale_kw=Dict("linthresh"=>0.01),xlocator=[0],xticklabels=["0"],
		yloc="neither",xloc="bottom",xlabel=L"$w_{wtg}$"
	)

	axs[9].plot(rh_rce,p_rce/100,c="gray3",linestyle="-.")
	axs[9].plot(rh_wet,p_dgw/100)
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
	fig.savefig(projectdir("figures","fig1-wtgschematic.png"),transparent=false,dpi=600)
	load(projectdir("figures","fig1-wtgschematic.png"))
end

# ╔═╡ Cell order:
# ╟─d41d44c4-e9d0-11ec-0b8d-0d5edc09dbd2
# ╟─b922b1af-31e6-4ea1-b30b-c5f03f882857
# ╟─181a3e2a-126b-4392-a2de-b2dbd68c9d15
# ╟─a29eda09-ac40-4131-8d86-e2ba55b50976
# ╠═a64e260a-9dd0-491e-b0fc-c4b6fdfa8b6f
# ╠═4e4bf212-220d-419e-9e90-fbd524c17881
# ╟─5eda6b78-22dc-4554-9a5d-ac213a9fa172
# ╟─0a73df1e-4244-47a0-bafb-b79c0a2d3357
# ╟─ae968145-c8b4-4e7c-b348-7d20085fac96
# ╟─c76f9180-d477-4ec4-84f5-ac2d6adb95b3
# ╟─16fb08d4-7bbc-4422-8755-530528e7a8f1
# ╟─41c12fde-938b-40f3-acf6-fecbae34e21c
