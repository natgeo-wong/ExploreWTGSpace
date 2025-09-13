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
# 02. WTG Adjustment Strengths
"

# ╔═╡ 3c16bd2b-55e0-492f-baf3-aef6a998e9d6
begin
	configDGW = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	configWTG = [
		0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),1,
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=5,nrows=2,aspect=0.4,axwidth=0.95,sharey=0,wspace=[1,2,1,1])

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conDGW in configDGW

		fnc = joinpath("DGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[8].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[8].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("DGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[3].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[3].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("KGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[9].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[9].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("KGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
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

		fnc = joinpath("TDG","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[10].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[10].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("TDG","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[5].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[5].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

	end
	
	for conWTG in configWTG

		fnc = joinpath("SPC","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				axs[7].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[7].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("SPC","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				axs[2].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[2].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("TGR","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				axs[6].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[6].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("TGR","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
	axs[2].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[3].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[4].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[5].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[6].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[7].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[8].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[9].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[10].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")

	axs[1].format(ylabel=L"$\tau$ / hr",ltitle="WTG")
	axs[2].format(ltitle="SWTG")
	axs[3].format(ltitle="DGW (Blossey)")
	axs[4].format(ltitle="DGW (Kuang)")
	axs[5].format(ylabel=L"$\alpha$",ltitle=L"DGW ($\partial_t\neq0$)")
	axs[6].format(ylabel=L"$\tau$ / hr")
	axs[10].format(ylabel=L"$\alpha$")

	for ax in axs
		ax.format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,10),xlabel=L"Hourly-Averaged Precipitation Rate / mm hr$^{-1}$",
			lrtitle="Wet",lltitle="Dry",urtitle="Wet",ultitle="Dry",yscale="log"
		)
	end

	for ii = vcat(1:2,6:7)
		axs[ii].format(ylim=(0.05,200))
	end
	for ii = vcat(3:5,8:10)
		axs[ii].format(ylim=(0.004,2500),ytickloc="r")
	end
	for ii = vcat(2:4,7:9)
		axs[ii].format(yticklabels=["","","",""],ytickminor=10:10:100,)
	end
	
	fig.savefig(plotsdir("02-wtgstrength.png"),transparent=false,dpi=400)
	load(plotsdir("02-wtgstrength.png"))
end
  ╠═╡ =#

# ╔═╡ 946c2a8f-6aa2-4fe9-a550-e2f84f1bdaee
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=4,nrows=2,aspect=0.4,axwidth=1,sharey=0,wspace=[1,2,1])

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conDGW in configDGW

		fnc = joinpath("DGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[7].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[7].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("DGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[3].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[3].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("KGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
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
				axs[8].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				axs[8].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("KGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
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

		fnc = joinpath("SPC","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				axs[6].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[6].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("SPC","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				axs[2].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[2].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("TGR","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				axs[5].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				axs[5].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("TGR","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		prcp  = ds_wtgprcp["precipitation"][:,:] / 24
		for ien = 1 : 15

			prcpii  = prcp[end-2399:end,ien]
			prcpii  = prcpii[.!isnan.(prcpii)]
			prcpμ   = mean(prcpii)
			prcpσ   = zeros(2,1)

			if !isnan(prcpμ)
				prcpσ[1] = abs(prcpμ - quantile(prcpii,0.05))
				prcpσ[2] = abs(quantile(prcpii,0.95) - prcpμ)
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
	axs[2].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[3].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[4].plot([1,1]*prcpRCEμ_P,[0,5000],c="grey")
	axs[5].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[6].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[7].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")
	axs[8].plot([1,1]*prcpRCEμ_T,[0,5000],c="grey")

	axs[1].format(ylabel=L"$\tau$ / hr",ltitle="(a) WTG-RZ")
	axs[2].format(ltitle="(b) SWTG",leftlabels=["RRTM","Idealized Radiative Cooling"])
	axs[3].format(ltitle="(c) DGW-B")
	axs[4].format(ltitle=L"(d) DGW-K ($\partial_t=0$)",ylabel=L"$\alpha$")
	axs[5].format(ylabel=L"$\tau$ / hr")
	axs[8].format(ylabel=L"$\alpha$")

	for ax in axs
		ax.format(
			xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
			xlim=(0,10),xlabel=L"Hourly-Averaged Precipitation Rate / mm hr$^{-1}$",
			lrtitle="Wet",lltitle="Dry",urtitle="Wet",ultitle="Dry",yscale="log"
		)
	end

	for ii = vcat(1:2,5:6)
		axs[ii].format(ylim=(0.05,200))
	end
	for ii = vcat(3:4,7:8)
		axs[ii].format(ylim=(0.004,2500),ytickloc="r")
	end
	for ii = vcat(2:3,6:7)
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
# ╟─57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╟─946c2a8f-6aa2-4fe9-a550-e2f84f1bdaee
