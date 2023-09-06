### A Pluto.jl notebook ###
# v0.19.27

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
	@quickactivate "2023GL104350"
	using NCDatasets
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the 2023GL104350 paper submission ..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# Figure S1. Vertical Velocity Profiles
"

# ╔═╡ 3c16bd2b-55e0-492f-baf3-aef6a998e9d6
begin
	configDGW1 = [0.2,0.5,1,2,5,10,20,50,100,200,500]
	configDGW2 = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	configWTG1 = [
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
	configWTG2 = [
		0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),1,
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
	nDGW1 = length(configDGW1)
	wet_DGW1 = pplt.get_colors("Blues",(nDGW1+2))
	dry_DGW1 = pplt.get_colors("Brown",(nDGW1+2))
	rce_DGW1 = pplt.get_colors("Teal", (nDGW1+2))
	nDGW2 = length(configDGW2)
	wet_DGW2 = pplt.get_colors("Blues",(nDGW2+2))
	dry_DGW2 = pplt.get_colors("Brown",(nDGW2+2))
	rce_DGW2 = pplt.get_colors("Teal", (nDGW2+2))
	nWTG1 = length(configWTG1)
	wet_WTG1 = pplt.get_colors("Blues",(nWTG1+2))
	dry_WTG1 = pplt.get_colors("Brown",(nWTG1+2))
	rce_WTG1 = pplt.get_colors("Teal", (nWTG1+2))
	nWTG2 = length(configWTG2)
	wet_WTG2 = pplt.get_colors("Blues",(nWTG2+2))
	dry_WTG2 = pplt.get_colors("Brown",(nWTG2+2))
	rce_WTG2 = pplt.get_colors("Teal", (nWTG2+2))
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=6,aspect=0.3,axwidth=1.1,wspace=[1,2.5,1,2.5,1])

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conii in 1 : nDGW1

		fnc = "DGW-P1282km300V64-$(dampingstrprnt(configDGW1[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))

		p    = ds_wwtg["p"][:]
		wwtg = ds_wwtg["wwtg"][:]*100
		prcp = ds_prcp["precipitation"][:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[1].plot(wwtgii,p[:,ien],c=wet_DGW1[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[1].plot(wwtgii,p[:,ien],c=dry_DGW1[conii+1])
			else
				axs[1].plot(wwtgii,p[:,ien],c=rce_DGW1[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

	end

	for conii in 1 : nDGW2

		fnc = "DGW-T1282km300V64-$(dampingstrprnt(configDGW2[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))

		p    = ds_wwtg["p"][:]
		wwtg = ds_wwtg["wwtg"][:]*100
		prcp = ds_prcp["precipitation"][:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[2].plot(wwtgii,p[:,ien],c=wet_DGW2[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[2].plot(wwtgii,p[:,ien],c=dry_DGW2[conii+1])
			else
				axs[2].plot(wwtgii,p[:,ien],c=rce_DGW2[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

	end

	for conii in 1 : nWTG1

		fnc = "TGR-P1282km300V64-$(relaxscalestrprnt(configWTG1[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))

		p    = ds_wwtg["p"][:]
		wwtg = ds_wwtg["wwtg"][:]*100
		prcp = ds_prcp["precipitation"][:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[3].plot(wwtgii,p[:,ien],c=wet_WTG1[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[3].plot(wwtgii,p[:,ien],c=dry_WTG1[conii+1])
			else
				axs[3].plot(wwtgii,p[:,ien],c=rce_WTG1[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		fnc = "SPC-P1282km300V64-$(relaxscalestrprnt(configWTG1[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))

		p    = ds_wwtg["p"][:]
		wwtg = ds_wwtg["wwtg"][:]*100
		prcp = ds_prcp["precipitation"][:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[5].plot(wwtgii,p[:,ien],c=wet_WTG1[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[5].plot(wwtgii,p[:,ien],c=dry_WTG1[conii+1])
			else
				axs[5].plot(wwtgii,p[:,ien],c=rce_WTG1[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

	end

	for conii in 1 : nWTG2

		fnc = "TGR-T1282km300V64-$(relaxscalestrprnt(configWTG2[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))

		p    = ds_wwtg["p"][:]
		wwtg = ds_wwtg["wwtg"][:]*100
		prcp = ds_prcp["precipitation"][:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[4].plot(wwtgii,p[:,ien],c=wet_WTG2[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[4].plot(wwtgii,p[:,ien],c=dry_WTG2[conii+1])
			else
				axs[4].plot(wwtgii,p[:,ien],c=rce_WTG2[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		fnc = "SPC-T1282km300V64-$(relaxscalestrprnt(configWTG2[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))

		p    = ds_wwtg["p"][:]
		wwtg = ds_wwtg["wwtg"][:]*100
		prcp = ds_prcp["precipitation"][:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[6].plot(wwtgii,p[:,ien],c=wet_WTG2[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[6].plot(wwtgii,p[:,ien],c=dry_WTG2[conii+1])
			else
				axs[6].plot(wwtgii,p[:,ien],c=rce_WTG2[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

	end

	axs[1].format(ltitle="(a) DGW")
	axs[3].format(ltitle="(b) TGR",)
	axs[5].format(ltitle="(c) SPC")

	for ax in axs
		ax.plot([0,0],[1000,10],c="k")
		ax.format(
			ylim=(1000,20),yscale="log",
			xscale="symlog",xscale_kw=Dict("linthresh"=>1),
			xlim=(-25,25),xlabel=L"$w_{wtg}$ / 10$^{-2}$ m s$^{-1}$",
		)
	end

	for ii in 1 : 2 : 5
		axs[ii].format(ultitle="(i) RRTM")
		axs[ii+1].format(ultitle="(ii) Ideal")
	end
	
	fig.savefig(projectdir("figures","figS1-wwtg.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS1-wwtg.png"))
end

# ╔═╡ Cell order:
# ╠═e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╟─57bbb1c8-e0f8-4b49-8f55-86af376c5167
