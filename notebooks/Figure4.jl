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
	
md"Loading modules for the 2023GL104350 project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# Figure 4. Baroclinic Mode Strengths
"

# ╔═╡ 0d1831f3-63e6-4467-b5cd-af86df6d8150
begin
    hvec = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,0.5,0.0]
	fvec = [0.0,0.1,0.2,0.3,0.4,0.5,1.0,1.0,1.0]
	nexp = length(hvec)
	md"List of baroclinic mode strength experiments"
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(
		[1,1,1,2,3,3,3,4],aspect=1.5,axwidth=2,wspace=[0,0,0,3,0,0,0],sharex=0
	)

	ds_rce = NCDataset(datadir("precipitation","RCE-T1282km300V64.nc"))
	prcp = mean(ds_rce["precipitation"][:]) / 24
	close(ds_rce)
	axs[1].plot([0,1],[1,1]*prcp,c="grey",lw=1)
	axs[2].plot([0,1],[1,1]*prcp,c="grey",lw=1)
	axs[3].plot([0,1],[1,1]*prcp,c="grey",lw=1)
	axs[4].plot([0,1],[1,1]*prcp,c="grey",lw=1)

	for iexp in 1 : (nexp-2)
		expname = "H$(@sprintf("%4.2f",hvec[iexp]))F$(@sprintf("%4.2f",fvec[iexp]))"
        expname = replace(expname,"."=>"d")

		ds_dgw  = NCDataset(datadir("precipitation","VDD-T1282km300V64-$expname.nc"))
		prcpdgw = ds_dgw["precipitation"][:]/24
		close(ds_dgw)

		ds_wtg  = NCDataset(datadir("precipitation","VDT-T1282km300V64-$expname.nc"))
		prcpwtg = ds_wtg["precipitation"][:]/24
		close(ds_wtg)
		
		prcpσ   = zeros(2,1)

		for ien = 1 : 15

			prcpii = prcpdgw[end-1000:end,ien]
			prcpμ  = mean(prcpii[.!isnan.(prcpii)])

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcp * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcp / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[1].scatter(fvec[iexp],prcpμ,c=mclr,s=20,zorder=5)
				axs[1].errorbar(fvec[iexp],prcpμ,prcpσ,0,c=lclr,zorder=4)
			end
			
			prcpii = prcpwtg[end-1000:end,ien]
			prcpμ  = mean(prcpii[.!isnan.(prcpii)])

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcp * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcp / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[3].scatter(fvec[iexp],prcpμ,c=mclr,s=20,zorder=5)
				axs[3].errorbar(fvec[iexp],prcpμ,prcpσ,0,c=lclr,zorder=4)
			end

		end
		
	end

	for iexp in (nexp-2)
		expname = "H$(@sprintf("%4.2f",hvec[iexp]))F$(@sprintf("%4.2f",fvec[iexp]))"
        expname = replace(expname,"."=>"d")

		ds_dgw  = NCDataset(datadir("precipitation","VDD-T1282km300V64-$expname.nc"))
		prcpdgw = ds_dgw["precipitation"][:]/24
		close(ds_dgw)

		ds_wtg  = NCDataset(datadir("precipitation","VDT-T1282km300V64-$expname.nc"))
		prcpwtg = ds_wtg["precipitation"][:]/24
		close(ds_wtg)
		
		prcpσ   = zeros(2,1)

		for ien = 1 : 15

			prcpii = prcpdgw[end-1000:end,ien]
			prcpμ  = mean(prcpii[.!isnan.(prcpii)])

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcp * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcp / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[1].scatter(0.6,prcpμ,c=mclr,s=20,zorder=5)
				axs[1].errorbar(0.6,prcpμ,prcpσ,0,c=lclr,zorder=4)
			end

			prcpii = prcpwtg[end-1000:end,ien]
			prcpμ  = mean(prcpii[.!isnan.(prcpii)])

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcp * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcp / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[3].scatter(0.6,prcpμ,c=mclr,s=20,zorder=5)
				axs[3].errorbar(0.6,prcpμ,prcpσ,0,c=lclr,zorder=4)
			end

		end
		
	end

	for iexp in (nexp-2) : nexp
		expname = "H$(@sprintf("%4.2f",hvec[iexp]))F$(@sprintf("%4.2f",fvec[iexp]))"
        expname = replace(expname,"."=>"d")

		ds_dgw  = NCDataset(datadir("precipitation","VDD-T1282km300V64-$expname.nc"))
		prcpdgw = ds_dgw["precipitation"][:]/24
		close(ds_dgw)

		ds_wtg  = NCDataset(datadir("precipitation","VDT-T1282km300V64-$expname.nc"))
		prcpwtg = ds_wtg["precipitation"][:]/24
		close(ds_wtg)
		
		prcpσ   = zeros(2,1)

		for ien = 1 : 15

			prcpii = prcpdgw[end-1000:end,ien]
			prcpμ  = mean(prcpii[.!isnan.(prcpii)])

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcp * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcp / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[2].scatter(hvec[iexp],prcpμ,c=mclr,s=20,zorder=5)
				axs[2].errorbar(hvec[iexp],prcpμ,prcpσ,0,c=lclr,zorder=4)
			end

			prcpii = prcpwtg[end-1000:end,ien]
			prcpμ  = mean(prcpii[.!isnan.(prcpii)])

			if !isnan(prcpμ)
				prcpσ[1] = prcpμ - quantile(prcpii,0.05)
				prcpσ[2] = quantile(prcpii,0.95) - prcpμ
				if prcpμ < prcp * 0.95
					mclr = "yellow9"
					lclr = "yellow3"
				elseif prcpμ > prcp / 0.95
					mclr = "blue9"
					lclr = "blue3"
				else
					mclr = "green9"
					lclr = "green3"
				end
				axs[4].scatter(hvec[iexp],prcpμ,c=mclr,s=20,zorder=5)
				axs[4].errorbar(hvec[iexp],prcpμ,prcpσ,0,c=lclr,zorder=4)
			end

		end
		
	end

	for iats = 1 : 2 : 4
		axs[iats].format(
			xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
		)
		axs[iats+1].format(
			yloc="right",xlim=(1,0),
			xlabel=L"\tau/\tau_h",xlocator=0:0.5:1,
		)
	end

	axs[1].format(
		xlabel=L"c_2",title=L"$a_m=10$ day$^{-1}$",
		ultitle=L"(a) $c_1 = 1$",
		xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
	)
	axs[2].format(
		xlim=(1,0),urtitle=L"(b) $c_2 = 1$",
		xlabel=L"c_1",xlocator=0:0.5:1,
	)

	axs[3].format(
		xlabel=L"c_2",title=L"$\tau=10$ hr",
		ultitle=L"(a) $c_1 = 1$",
		xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
	)
	axs[4].format(
		xlim=(1,0),urtitle=L"(b) $c_2 = 1$",
		xlabel=L"c_1",xlocator=0:0.5:1,
	)
	axs[2].format(ylim=(0,1),ylabel=L"Rain Rate / mm hr$^{-1}$")
	
	fig.savefig(projectdir("figures","fig4-baroclinicmodes.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig4-baroclinicmodes.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─0d1831f3-63e6-4467-b5cd-af86df6d8150
# ╟─55230f4a-7661-11eb-1c37-8b022b95e08e
