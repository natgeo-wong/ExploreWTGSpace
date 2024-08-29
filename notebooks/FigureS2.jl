### A Pluto.jl notebook ###
# v0.19.42

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
# Figure S2. Time Series, Full Radiation
"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configDGW = [0.02,0.05,0.1,0.2,500]
	configWTG = [0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),50*sqrt(2)]
	nconDGW = length(configDGW)
	nconWTG = length(configWTG)
	blues_DGW = pplt.get_colors("Blues",(nconDGW+2))
	blues_WTG = pplt.get_colors("Blues",(nconWTG+2))
	lgd_DGW = Dict("frame"=>false,"ncols"=>2)
	lgd_WTG = Dict("frame"=>false,"ncols"=>2)
	md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 864f1d31-c629-412b-825a-98fe9398b591
begin
	pplt.close()
	fts,ats = pplt.subplots(nrows=2,aspect=3,axwidth=4)

	for ic in 1 : nconDGW

		fnc = "DGW-D1282km300V64-$(dampingstrprnt(configDGW[ic])).nc"
		ds_dgwprcp = NCDataset(datadir("precipitation",fnc))

		tdgw    = ds_dgwprcp["time"][:]; #tdgw = reshape(tdgw,24,:)
		# tdgw    = dropdims(mean(tdgw,dims=1),dims=1)
		prcpdgw = ds_dgwprcp["precipitation"][:,:] / 24
		
		
		for ien = 1 : 15

			prcpii = prcpdgw[:,ien]
			# prcpii = reshape(prcpii,24,:)
			# prcpii = dropdims(mean(prcpii,dims=1),dims=1)

			if ien == 1
				constr = @sprintf("%.1e",configDGW[ic])
				ats[1].plot(
					tdgw,prcpii,color=blues_DGW[ic+1],
					label=(L"$\alpha =$" * " $(constr)"),
					legend="r",legend_kw=lgd_DGW
				)
			else
				ats[1].plot(tdgw,prcpii,color=blues_DGW[ic+1])
			end
			
		end

		close(ds_dgwprcp)

	end

	for ic in 1 : nconWTG

		# fnc = "TGR-P1282km300V64-$(relaxscalestrprnt(configWTG[ic])).nc"
		# ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		# twtg    = ds_wtgprcp["time"][:]; twtg = reshape(twtg,24,:)
		# twtg    = dropdims(mean(twtg,dims=1),dims=1)
		# prcpwtg = ds_wtgprcp["precipitation"][:] / 24
		
		
		# for ien = 1 : 15

		# 	prcpii = prcpwtg[:,ien]
		# 	prcpii = reshape(prcpii,24,:)
		# 	prcpii = dropdims(mean(prcpii,dims=1),dims=1)

		# 	if ien == 1
		# 		constr = @sprintf("%.1e",configWTG[ic])
		# 		ats[2].plot(
		# 			twtg,prcpii,color=blues_WTG[ic+1],
		# 			label=(L"$\tau =$" * " $(constr) hr"),
		# 			legend="r",legend_kw=lgd_WTG
		# 		)
		# 	else
		# 		ats[2].plot(twtg,prcpii,color=blues_WTG[ic+1])
		# 	end
			
		# end

		# close(ds_wtgprcp)


		fnc = "SPC-D1282km300V64-$(relaxscalestrprnt(configWTG[ic])).nc"
		ds_wtgprcp = NCDataset(datadir("precipitation",fnc))

		twtg    = ds_wtgprcp["time"][:]; #twtg = reshape(twtg,24,:)
		# twtg    = dropdims(mean(twtg,dims=1),dims=1)
		prcpwtg = ds_wtgprcp["precipitation"][:,:] / 24


		for ien = 1 : 1

			prcpii = prcpwtg[:,ien]
			# prcpii = reshape(prcpii,24,:)
			# prcpii = dropdims(mean(prcpii,dims=1),dims=1)

			if ien == 1
				constr = @sprintf("%.1e",configWTG[ic])
				ats[2].plot(
					twtg,prcpii,color=blues_WTG[ic+1],
					label=(L"$\tau =$" * " $(constr) hr"),
					legend="r",legend_kw=lgd_WTG
				)
			else
				ats[2].plot(twtg,prcpii,color=blues_WTG[ic+1])
			end

		end

		close(ds_wtgprcp)

	end

	ds_rceprcp = NCDataset(datadir("precipitation","RCE-P1282km300V64.nc"))

	# t_RCE = ds_rceprcp["time"][:]
	# prcp_RCE = ds_rceprcp["precipitation"][:] / 24
	# ats[1].plot(t_RCE,prcp_RCE[:,1],c="k",label="RCE",legend="r")
	# ats[2].plot(t_RCE,prcp_RCE[:,1],c="k",label="RCE",legend="r")
	# ats[3].plot(t_RCE,prcp_RCE[:,1],c="k",label="RCE",legend="r")
	# ats[1].plot(t_RCE,prcp_RCE[:,2:10],c="k")
	# ats[2].plot(t_RCE,prcp_RCE[:,2:10],c="k")
	# ats[3].plot(t_RCE,prcp_RCE[:,2:10],c="k")

	ats[1].format(ultitle="(a) Damped Gravity Wave")
	ats[2].format(ultitle="(b) Temperature Gradient Relaxation")
	ats[2].format(ultitle="(c) Spectral WTG")

	close(ds_rceprcp)
	
	for ax in ats
		ax.format(
			xlim=(0,250),#yscale="symlog",yscale_kw=Dict("linthresh"=>0.1),
			ylim=(0,5),
			ylabel=L"Daily-Averaged Rainfall Rate / mm hr$^{-1}$",xlabel="Days"
		)
	end
	
	fts.savefig(projectdir("figures","figS2-timeseries-P1282km300V64.png"),transparent=false,dpi=400)
	load(projectdir("figures","figS2-timeseries-P1282km300V64.png"))
end

# ╔═╡ 92277ba1-1cd2-4f8a-8587-1f32f04f068b


# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╠═864f1d31-c629-412b-825a-98fe9398b591
# ╠═92277ba1-1cd2-4f8a-8587-1f32f04f068b
