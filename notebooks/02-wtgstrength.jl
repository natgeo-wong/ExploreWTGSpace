### A Pluto.jl notebook ###
# v0.19.41

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
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace paper submission ..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 02. WTG Adjustment Strengths
"

# ╔═╡ 3c16bd2b-55e0-492f-baf3-aef6a998e9d6
begin
	configDGW = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	configDG2 = [0.02,0.05,0.1,0.2,0.5,1,2,5]
	configWTG = [
		0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),1,
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=6,aspect=0.3,axwidth=0.9,sharey=0,wspace=[1,1,4,1,1])

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-D1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_D = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conDGW in configDGW

		fnc = "DGW-T1282km300V64-$(dampingstrprnt(conDGW)).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_dgwprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcpRCEμ_T * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcpRCEμ_T / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[6].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[6].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = "DGW-D1282km300V64-$(dampingstrprnt(conDGW)).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_dgwprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcpRCEμ_D * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcpRCEμ_D / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[5].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[5].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = "DGW-P1282km300V64-$(dampingstrprnt(conDGW)).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_dgwprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcpRCEμ_P * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcpRCEμ_P / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[4].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[4].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

	end
	
	for conWTG in configWTG

		fnc = "SPC-T1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcpRCEμ_T * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcpRCEμ_T / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[3].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[3].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		fnc = "SPC-D1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcpRCEμ_D * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcpRCEμ_D / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[2].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[2].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = "SPC-P1282km300V64-$(relaxscalestrprnt(conWTG)).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcpRCEμ_P * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcpRCEμ_P / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[1].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[1].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

	end
	
	axs[1].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[2].plot([1,1]*prcpRCEμ_D,[0,5000],c="grey")
	axs[3].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[4].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[5].plot([1,1]*prcpRCEμ_D,[0,5000],c="grey")
	axs[6].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")

	axs[1].format(ultitle="Perpetual",ylabel=L"$\tau$ / hr",ltitle="Spectral WTG")
	axs[2].format(ultitle="Diurnal")
	axs[3].format(ultitle="Idealized")
	axs[4].format(ultitle="Perpetual",ltitle="Damped Gravity Wave")
	axs[5].format(ultitle="Diurnal")
	axs[6].format(ultitle="Idealized",ylabel=L"$\alpha$")

	for ax in axs
		ax.format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,2),xlabel=L"Hourly-Averaged Precipitation Rate / mm hr$^{-1}$",
			lrtitle="Wet",lltitle="Dry",yscale="log"
		)
	end

	for ii = 1 : 3
		axs[ii].format(ylim=(0.05,200))
	end
	for ii = 4 : 6
		axs[ii].format(ylim=(0.004,2500),ytickloc="r")
	end
	for ii = 2 : 5
		axs[ii].format(yticklabels=["","","",""],ytickminor=10:10:100,)
	end
	
	fig.savefig(plotsdir("02-wtgstrength.png"),transparent=false,dpi=400)
	load(plotsdir("02-wtgstrength.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╠═57bbb1c8-e0f8-4b49-8f55-86af376c5167
