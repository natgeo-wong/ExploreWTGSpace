### A Pluto.jl notebook ###
# v0.20.5

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
	using NCDatasets
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("ultraplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace paper submission ..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# Figure S1. Vertical Velocity Profiles
"

# ╔═╡ 3c16bd2b-55e0-492f-baf3-aef6a998e9d6
begin
	configDGW = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	configWTG = [
		0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),1,
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
	nDGW = length(configDGW)
	wet_DGW = pplt.get_colors("Blues",(nDGW+2))
	dry_DGW = pplt.get_colors("Brown",(nDGW+2))
	rce_DGW = pplt.get_colors("Teal", (nDGW+2))
	nWTG = length(configWTG)
	wet_WTG = pplt.get_colors("Blues",(nWTG+2))
	dry_WTG = pplt.get_colors("Brown",(nWTG+2))
	rce_WTG = pplt.get_colors("Teal", (nWTG+2))
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=4,nrows=2,aspect=0.4,axwidth=1,wspace=1,)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conii in 1 : nDGW

		fnc = "$(dampingstrprnt(configDGW[conii])).nc"
		ds_wwtg = NCDataset(datadir("wwtg","DGW","P1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","DGW","P1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[3].plot(wwtgii,p[:,ien],c=wet_DGW[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[3].plot(wwtgii,p[:,ien],c=dry_DGW[conii+1])
			else
				axs[3].plot(wwtgii,p[:,ien],c=rce_DGW[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)
		
		ds_wwtg = NCDataset(datadir("wwtg","KGW","P1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","KGW","P1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[4].plot(wwtgii,p[:,ien],c=wet_DGW[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[4].plot(wwtgii,p[:,ien],c=dry_DGW[conii+1])
			else
				axs[4].plot(wwtgii,p[:,ien],c=rce_DGW[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		ds_wwtg = NCDataset(datadir("wwtg","DGW","T1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","DGW","T1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[7].plot(wwtgii,p[:,ien],c=wet_DGW[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[7].plot(wwtgii,p[:,ien],c=dry_DGW[conii+1])
			else
				axs[7].plot(wwtgii,p[:,ien],c=rce_DGW[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		close(ds_wwtg)
		close(ds_prcp)

		ds_wwtg = NCDataset(datadir("wwtg","KGW","T1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","KGW","T1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[8].plot(wwtgii,p[:,ien],c=wet_DGW[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[8].plot(wwtgii,p[:,ien],c=dry_DGW[conii+1])
			else
				axs[8].plot(wwtgii,p[:,ien],c=rce_DGW[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

	end

	for conii in 1 : nWTG

		fnc = "$(relaxscalestrprnt(configWTG[conii])).nc"
		
		ds_wwtg = NCDataset(datadir("wwtg","TGR","P1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","TGR","P1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[1].plot(wwtgii,p[:,ien],c=wet_WTG[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[1].plot(wwtgii,p[:,ien],c=dry_WTG[conii+1])
			else
				axs[1].plot(wwtgii,p[:,ien],c=rce_WTG[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		ds_wwtg = NCDataset(datadir("wwtg","SPC","P1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","SPC","P1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.95
				axs[2].plot(wwtgii,p[:,ien],c=wet_WTG[conii+1])
			elseif prcpμ < prcpRCEμ_P * 0.95
				axs[2].plot(wwtgii,p[:,ien],c=dry_WTG[conii+1])
			else
				axs[2].plot(wwtgii,p[:,ien],c=rce_WTG[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		ds_wwtg = NCDataset(datadir("wwtg","TGR","T1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","TGR","T1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[5].plot(wwtgii,p[:,ien],c=wet_WTG[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[5].plot(wwtgii,p[:,ien],c=dry_WTG[conii+1])
			else
				axs[5].plot(wwtgii,p[:,ien],c=rce_WTG[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

		ds_wwtg = NCDataset(datadir("wwtg","SPC","T1282km300V64",fnc))
		ds_prcp = NCDataset(datadir("precipitation","SPC","T1282km300V64",fnc))

		p    = ds_wwtg["p"][:,:]
		wwtg = ds_wwtg["wwtg"][:,:,:]*100
		prcp = ds_prcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.95
				axs[6].plot(wwtgii,p[:,ien],c=wet_WTG[conii+1])
			elseif prcpμ < prcpRCEμ_T * 0.95
				axs[6].plot(wwtgii,p[:,ien],c=dry_WTG[conii+1])
			else
				axs[6].plot(wwtgii,p[:,ien],c=rce_WTG[conii+1])
			end
			
		end

		close(ds_wwtg)
		close(ds_prcp)

	end

	axs[1].format(ltitle="(a) WTG-RZ")
	axs[2].format(ltitle="(b) SWTG",)
	axs[3].format(ltitle="(c) DGW-B")
	axs[4].format(ltitle=L"(d) DGW-K ($\partial_t = 0$)")

	for ax in axs
		ax.plot([0,0],[1000,10],c="k")
		ax.format(
			ylim=(1000,20),yscale="log",ylabel="Pressure / hPa",
			xscale="symlog",xscale_kw=Dict("linthresh"=>1),
			xlim=(-25,25),xlabel=L"$w_{wtg}$ / 10$^{-2}$ m s$^{-1}$",
			rightlabels=["RRTM","Idealized Radiative Cooling"]
		)
	end
	
	fig.savefig(projectdir("figures","figS1-wwtg.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS1-wwtg.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╟─57bbb1c8-e0f8-4b49-8f55-86af376c5167
