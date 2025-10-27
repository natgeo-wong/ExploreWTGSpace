### A Pluto.jl notebook ###
# v0.20.19

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
	using Printf
	using NCDatasets, StatsBase
	using CairoMakie, MakieExtra, LaTeXStrings
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
	fig = Figure()
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 87038073-d846-40c2-8a0f-b5969185998c
function scattercolor(μ,RCE)

	if μ < RCE * 0.95
		mclr = :gold2
		lclr = :gold2
	elseif μ > RCE / 0.95
		mclr = :dodgerblue2
		lclr = :dodgerblue2
	else
		mclr = :palegreen3
		lclr = :palegreen3
	end

	return mclr,lclr

end

# ╔═╡ 946c2a8f-6aa2-4fe9-a550-e2f84f1bdaee
begin
	empty!(fig)
	
	axs = Matrix{Axis}(undef,2,4)

	for icol = 1 : 4, irow = 1 : 2
		if isone(irow)
			axs[irow,icol] = Axis(
				fig[irow,icol],width=100,height=250,
				xscale=SymLog(0.01),xticks=[0,0.01,0.1,1,10],
				yticklabelsvisible=isone(icol)||icol==4,
				yscale=log10,yticks=[0.01,0.1,1,10,100,1000],
				xminorticks=vcat((1:9)/100,(1:9)/10,1:9),
				xminorticksvisible=true,xticklabelsvisible=false
			)
		else
			axs[irow,icol] = Axis(
				fig[irow,icol],width=100,height=250,
				xscale=SymLog(0.01),xticks=[0,0.01,0.1,1,10],
				yticklabelsvisible=isone(icol)||icol==4,
				yscale=log10,yticks=[0.01,0.1,1,10,100,1000],
				xminorticks=vcat((1:9)/100,(1:9)/10,1:9),
				xminorticksvisible=true
			)
		end
	end

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","D0642km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","T0642km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:] / 24
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conDGW in configDGW

		fnc = joinpath("DGW","T0642km300V64","$(dampingstrprnt(conDGW)).nc")
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_T)
				scatter!(axs[2,3],prcpμ,conDGW,color=mclr,markersize=12)
				# axs[7].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

		fnc = joinpath("DGW","D0642km300V64","$(dampingstrprnt(conDGW)).nc")
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_P)
				scatter!(axs[1,3],prcpμ,conDGW,color=mclr,markersize=12)
				# axs[3].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				# axs[3].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_T)
				scatter!(axs[2,4],prcpμ,conDGW,color=mclr,markersize=12)
				# axs[8].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				# axs[8].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_P)
				scatter!(axs[1,4],prcpμ,conDGW,color=mclr,markersize=12)
				# axs[4].scatter(prcpμ,conDGW,c=mclr,s=20,zorder=5)
				# axs[4].errorbar(prcpμ,conDGW,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_dgwprcp)

	end
	
	for conWTG in configWTG

		fnc = joinpath("SPC","T0642km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_T)
				scatter!(axs[2,2],prcpμ,conWTG,color=mclr,markersize=12)
				# axs[6].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				# axs[6].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

		fnc = joinpath("SPC","D0642km300V64","$(relaxscalestrprnt(conWTG)).nc")
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_P)
				scatter!(axs[1,2],prcpμ,conWTG,color=mclr,markersize=12)
				# axs[2].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				# axs[2].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_T)
				scatter!(axs[2,1],prcpμ,conWTG,color=mclr,markersize=12)
				# axs[5].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				# axs[5].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
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
				mclr,lclr = scattercolor(prcpμ,prcpRCEμ_P)
				scatter!(axs[1,1],prcpμ,conWTG,color=mclr,markersize=12)
				# axs[1].scatter(prcpμ,conWTG,c=mclr,s=20,zorder=5)
				# axs[1].errorbar(prcpμ,conWTG,0,prcpσ,c=lclr,zorder=4)
			end
			
		end

		close(ds_wtgprcp)

	end
	
	lines!(axs[1,1],[1,1]*prcpRCEμ_P,[0.002,5000],color=:grey)
	lines!(axs[1,2],[1,1]*prcpRCEμ_P,[0.002,5000],color=:grey)
	lines!(axs[1,3],[1,1]*prcpRCEμ_P,[0.002,5000],color=:grey)
	lines!(axs[1,4],[1,1]*prcpRCEμ_P,[0.002,5000],color=:grey)
	lines!(axs[2,1],[1,1]*prcpRCEμ_T,[0.002,5000],color=:grey)
	lines!(axs[2,2],[1,1]*prcpRCEμ_T,[0.002,5000],color=:grey)
	lines!(axs[2,3],[1,1]*prcpRCEμ_T,[0.002,5000],color=:grey)
	lines!(axs[2,4],[1,1]*prcpRCEμ_T,[0.002,5000],color=:grey)

	for icol = 1 : 4, irow = 1 : 2
		xlims!(axs[irow,icol],0,5)
		if icol >= 3
			ylims!(axs[irow,icol],0.004,2500)
			axs[irow,icol].yaxisposition=:right
		else
			ylims!(axs[irow,icol],0.05,200)
		end
	end

	# axs[1].format(ylabel=L"$\tau$ / hr",ltitle="(a) WTG-RZ")
	# axs[2].format(ltitle="(b) SWTG",leftlabels=["RRTM","Idealized Radiative Cooling"])
	# axs[3].format(ltitle="(c) DGW-B")
	# axs[4].format(ltitle=L"(d) DGW-K ($\partial_t=0$)",ylabel=L"$\alpha$")
	# axs[5].format(ylabel=L"$\tau$ / hr")
	# axs[8].format(ylabel=L"$\alpha$")

	# for ax in axs
	# 	ax.format(
	# 		xscale="symlog",xscale_kw=Dict("linthresh"=>0.01),
	# 		xlim=(0,10),xlabel=L"Hourly-Averaged Precipitation Rate / mm hr$^{-1}$",
	# 		lrtitle="Wet",lltitle="Dry",urtitle="Wet",ultitle="Dry",yscale="log"
	# 	)
	# end

	# for ii = vcat(1:2,5:6)
	# 	axs[ii].format(ylim=(0.05,200))
	# end
	# for ii = vcat(3:4,7:8)
	# 	axs[ii].format(ylim=(0.004,2500),ytickloc="r")
	# end
	# for ii = vcat(2:3,6:7)
	# 	axs[ii].format(yticklabels=["","","",""],ytickminor=10:10:100,)
	# end

	Label(fig[3,:],L"Precipitation Rate / mm hr$^{-1}$")
	Label(fig[:,0],L"$\tau$ / hr",rotation=π/2)
	Label(fig[:,5],L"$a_m$ / day$^{-1}$, $\lambda/4 =$ 650 km",rotation=π/2)

	resize_to_layout!(fig)
	fig
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╟─87038073-d846-40c2-8a0f-b5969185998c
# ╟─946c2a8f-6aa2-4fe9-a550-e2f84f1bdaee
