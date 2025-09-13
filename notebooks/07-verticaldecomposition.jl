### A Pluto.jl notebook ###
# v0.19.40

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
	nDGW2 = length(configDGW2)
	nWTG1 = length(configWTG1)
	nWTG2 = length(configWTG2)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 8cd96b4d-58f5-4cab-8fb0-8fe4a6cd5f09
function vmcoef(w,z;ncoef=20)

	coefmat = zeros(ncoef)

	ivec = .!iszero.(w); nz = sum(ivec)+1; ztrop = z[nz]
	wvec = vcat(0,w[ivec],0)
	zvec = vcat(0,z[ivec],ztrop)
	

	for icoef = 1 : ncoef, iz = 1 : nz

		iiz = (zvec[iz+1] + zvec[iz]) / 2
		iiw =  wvec[iz+1] + wvec[iz]
		coefmat[icoef] += iiw * sin(icoef*pi*iiz/ztrop) * (iiz-zvec[iz])

	end

	return coefmat / ztrop

end

# ╔═╡ aafb9ac1-561e-49d3-9136-d509ef17e5d5
begin
	pplt.close()
	ftmp,atmp = pplt.subplots(ncols=4,nrows=2,axwidth=1,wspace=[1,3,1])
	lvls = -0.5:0.05:0.5

	nmodes = 32
	wetmat1 = zeros(nDGW1,nmodes)
	drymat1 = zeros(nDGW1,nmodes)
	wetmat2 = zeros(nDGW2,nmodes)
	drymat2 = zeros(nDGW2,nmodes)
		
	ds_rceprcp = NCDataset(datadir("precipitation","RCE-P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:]
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)
		
	ds_rceprcp = NCDataset(datadir("precipitation","RCE-T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:]
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for iconDGW = 1 : nDGW1
		
		fnc = "DGW-P1282km300V64-$(dampingstrprnt(configDGW1[iconDGW])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))
		prcp = ds_prcp["precipitation"][:,:]
		z    = ds_wwtg["z"][:]
		wwtg = ds_wwtg["wwtg"][:,:,:]
		close(ds_wwtg)
		for ien = 1 : 15
			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_P / 0.9
				wetmat1[iconDGW,:] += vmcoef(wwtgii,z,ncoef=nmodes)
			elseif prcpμ < prcpRCEμ_P * 0.9
				drymat1[iconDGW,:] += vmcoef(wwtgii,z,ncoef=nmodes)
			end
		end
	end
	
	for iconDGW = 1 : nDGW2
		
		fnc = "DGW-T1282km300V64-$(dampingstrprnt(configDGW2[iconDGW])).nc"
		ds_wwtg = NCDataset(datadir("wwtg",fnc))
		ds_prcp = NCDataset(datadir("precipitation",fnc))
		prcp = ds_prcp["precipitation"][:,:]
		z    = ds_wwtg["z"][:]
		wwtg = ds_wwtg["wwtg"][:,:,:]
		close(ds_wwtg)
		for ien = 1 : 15
			wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
			prcpμ  = mean(prcp[end-2399:end,ien])
			if prcpμ > prcpRCEμ_T / 0.9
				wetmat2[iconDGW,:] += vmcoef(wwtgii,z,ncoef=nmodes)
			elseif prcpμ < prcpRCEμ_T * 0.9
				drymat2[iconDGW,:] += vmcoef(wwtgii,z,ncoef=nmodes)
			end
		end
	end
	
	drymat1 ./= abs.(drymat1[:,1])
	wetmat1 ./= abs.(wetmat1[:,1])
	drymat2 ./= abs.(drymat2[:,1])
	wetmat2 ./= abs.(wetmat2[:,1])
	
	c = 
	atmp[1].pcolormesh(1:nmodes,configDGW1,drymat1,levels=lvls,extend="both")
	atmp[2].pcolormesh(1:nmodes,configDGW2,drymat2,levels=lvls,extend="both")
	atmp[5].pcolormesh(1:nmodes,configDGW1,wetmat1,levels=lvls,extend="both")
	atmp[6].pcolormesh(1:nmodes,configDGW2,wetmat2,levels=lvls,extend="both")

	# wetmat1 = zeros(nWTG1,nmodes)
	# drymat1 = zeros(nWTG1,nmodes)
	# wetmat2 = zeros(nWTG2,nmodes)
	# drymat2 = zeros(nWTG2,nmodes)

	# for iconWTG = 1 : nWTG1
		
	# 	fnc = "SPC-P1282km300V64-$(relaxscalestrprnt(configWTG1[iconWTG])).nc"
	# 	ds_wwtg = NCDataset(datadir("wwtg",fnc))
	# 	ds_prcp = NCDataset(datadir("precipitation",fnc))
	# 	prcp = ds_prcp["precipitation"][:,:]
	# 	z    = ds_wwtg["z"][:]
	# 	wwtg = ds_wwtg["wwtg"][:,:,:]
	# 	close(ds_wwtg)
	# 	for ien = 1 : 15
	# 		wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
	# 		prcpμ  = mean(prcp[end-2399:end,ien])
	# 		if prcpμ > prcpRCEμ_P / 0.9
	# 			wetmat1[iconWTG,:] += vmcoef(wwtgii,z,ncoef=nmodes)
	# 		elseif prcpμ < prcpRCEμ_P * 0.9
	# 			drymat1[iconWTG,:] += vmcoef(wwtgii,z,ncoef=nmodes)
	# 		end
	# 	end
	# end
	
	# for iconWTG = 1 : nWTG2
		
	# 	fnc = "SPC-T1282km300V64-$(relaxscalestrprnt(configWTG2[iconWTG])).nc"
	# 	ds_wwtg = NCDataset(datadir("wwtg",fnc))
	# 	ds_prcp = NCDataset(datadir("precipitation",fnc))
	# 	prcp = ds_prcp["precipitation"][:,:]
	# 	z    = ds_wwtg["z"][:]
	# 	wwtg = ds_wwtg["wwtg"][:,:,:]
	# 	close(ds_wwtg)
	# 	for ien = 1 : 15
	# 		wwtgii = dropdims(mean(wwtg[:,end-2399:end,ien],dims=2),dims=2)
	# 		prcpμ  = mean(prcp[end-2399:end,ien])
	# 		if prcpμ > prcpRCEμ_T / 0.9
	# 			wetmat2[iconWTG,:] += vmcoef(wwtgii,z,ncoef=nmodes)
	# 		elseif prcpμ < prcpRCEμ_T * 0.9
	# 			drymat2[iconWTG,:] += vmcoef(wwtgii,z,ncoef=nmodes)
	# 		end
	# 	end
	# end
	
	# drymat1 ./= abs.(drymat1[:,1])
	# wetmat1 ./= abs.(wetmat1[:,1])
	# drymat2 ./= abs.(drymat2[:,1])
	# wetmat2 ./= abs.(wetmat2[:,1])

	# atmp[3].pcolormesh(1:nmodes,configWTG1,drymat1,levels=lvls,extend="both")
	# atmp[4].pcolormesh(1:nmodes,configWTG2,drymat2,levels=lvls,extend="both")
	# atmp[7].pcolormesh(1:nmodes,configWTG1,wetmat1,levels=lvls,extend="both")
	# atmp[8].pcolormesh(1:nmodes,configWTG2,wetmat2,levels=lvls,extend="both")

	for ax in atmp
		ax.format(yscale="log",ylim=(1,100),xlim=(0,30))
	end
	# atmp[1].format(yscale="log",ylim=(1,100),xlim=(0,6))
	# atmp[2].format(yscale="log",ylim=(1,100),xlim=(0,6))
	# atmp[5].format(yscale="log",ylim=(1,100),xlim=(0,6))
	# atmp[6].format(yscale="log",ylim=(1,100),xlim=(0,6))
	
	# atmp[3].format(yscale="log",ylim=(2,50),xlim=(0,6))
	# atmp[4].format(yscale="log",ylim=(2,50),xlim=(0,6))
	# atmp[7].format(yscale="log",ylim=(2,50),xlim=(0,6))
	# atmp[8].format(yscale="log",ylim=(2,50),xlim=(0,6))

	ftmp.colorbar(c)
	ftmp.savefig("test.png",transparent=false,dpi=200)
	load("test.png")
end

# ╔═╡ 88ac458d-f4a6-431e-ae69-1b8d70f3ee38
drymat1[:,1:6]

# ╔═╡ 57bbb1c8-e0f8-4b49-8f55-86af376c5167
# ╠═╡ show_logs = false
# ╠═╡ disabled = true
#=╠═╡
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
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─3c16bd2b-55e0-492f-baf3-aef6a998e9d6
# ╟─8cd96b4d-58f5-4cab-8fb0-8fe4a6cd5f09
# ╠═aafb9ac1-561e-49d3-9136-d509ef17e5d5
# ╠═88ac458d-f4a6-431e-ae69-1b8d70f3ee38
# ╟─57bbb1c8-e0f8-4b49-8f55-86af376c5167
