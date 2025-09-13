### A Pluto.jl notebook ###
# v0.20.0

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
	using MultivariateStats
	using NCDatasets
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace paper ..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 3. Time-Series and Power Spectrum

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ c0205377-e7c3-4272-8dc0-68230da2118b
begin
	configDGW = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
	configWTG = [
		0.1*sqrt(2),0.2,0.2*sqrt(2.5),0.5,0.5*sqrt(2),1,
		sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2),10,
		10*sqrt(2),20,20*sqrt(2.5),50,50*sqrt(2)
	]
	nmode = 10
	blues_WTG = pplt.get_colors("Blues_r",(nmode+2))
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 4a85fc54-279c-449d-b84b-1dcb00d1ca74
# ╠═╡ show_logs = false
begin
	pplt.close()
	fig,axs = pplt.subplots(ncols=5,nrows=2,aspect=0.4,axwidth=0.95,sharey=0,wspace=[1,2,1,1])

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","P1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:]
	prcpRCEμ_P = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	ds_rceprcp = NCDataset(datadir("precipitation","RCE","T1282km300V64.nc"))
	prcp_RCE = ds_rceprcp["precipitation"][:,:]
	prcpRCEμ_T = mean(prcp_RCE[1001:2000,:])
	close(ds_rceprcp)

	for conWTG in configWTG

		fnc = joinpath("TGR","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[1].scatter(-mean(wₙ[imode,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[1].scatter(mean(wₙ[imode,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("TGR","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[6].scatter(-mean(wₙ[imode,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[6].scatter(mean(wₙ[imode,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end

		fnc = joinpath("SPC","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[2].scatter(-mean(wₙ[imode,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[2].scatter(mean(wₙ[imode,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("SPC","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[7].scatter(-mean(wₙ[imode,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[7].scatter(mean(wₙ[imode,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end

	end

	for conDGW in configDGW

		fnc = joinpath("DGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[3].scatter(-mean(wₙ[imode,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[3].scatter(mean(wₙ[imode,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("DGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[8].scatter(-mean(wₙ[imode,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[8].scatter(mean(wₙ[imode,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end

		fnc = joinpath("KGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[4].scatter(-mean(wₙ[imode,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[4].scatter(mean(wₙ[imode,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("KGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[9].scatter(-mean(wₙ[imode,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[9].scatter(mean(wₙ[imode,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end

		fnc = joinpath("TDG","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[5].scatter(-mean(wₙ[imode,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[5].scatter(mean(wₙ[imode,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("TDG","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				axs[10].scatter(-mean(wₙ[imode,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				axs[10].scatter(mean(wₙ[imode,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end

	end

	axs[1].format(ylabel=L"$\tau$ / hr",ltitle="WTG")
	axs[2].format(ltitle="SWTG")
	axs[3].format(ltitle="DGW (Blossey)")
	axs[4].format(ltitle="DGW (Kuang)")
	axs[5].format(ylabel=L"$\alpha$",ltitle=L"DGW ($\partial_t\neq0$)")
	axs[6].format(ylabel=L"$\tau$ / hr")
	axs[10].format(ylabel=L"$\alpha$")

	for ax in axs
		ax.format(
			xlim=(-1,1),yscale="log",
			lrtitle="Wet",lltitle="Dry",urtitle="Wet",ultitle="Dry",
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
	
	fig.savefig(plotsdir("12-wn.png"),transparent=false,dpi=400)
	load(plotsdir("12-wn.png"))
end

# ╔═╡ f032247e-7959-4a7a-8de2-7e09cf918b48
# ╠═╡ show_logs = false
begin
	pplt.close()
	fnm,anm = pplt.subplots(ncols=5,nrows=2,aspect=0.4,axwidth=0.95,sharey=0,wspace=[1,2,1,1])

	for conWTG in configWTG

		fnc = joinpath("TGR","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[1].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[1].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("TGR","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[6].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[6].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end

		fnc = joinpath("SPC","P1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[2].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[2].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("SPC","T1282km300V64","$(relaxscalestrprnt(conWTG)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[7].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conWTG,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[7].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conWTG,c=blues_WTG[imode+1],s=10)
			end
		end

	end

	for conDGW in configDGW

		fnc = joinpath("DGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[3].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[3].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("DGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[8].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[8].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end

		fnc = joinpath("KGW","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[4].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[4].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("KGW","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[9].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[9].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end

		fnc = joinpath("TDG","P1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_P / 0.9)
		jj = prcp .> (prcpRCEμ_P * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[5].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[5].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end
		
		fnc = joinpath("TDG","T1282km300V64","$(dampingstrprnt(conDGW)).nc")
		
		ds   = NCDataset(datadir("precipitation",fnc))
		prcp = ds["precipitation"][end-2399:end,:]
		close(ds)
		ds   = NCDataset(datadir("wwtg",fnc))
		wₙ   = ds["wₙ"][:,:,:]
		close(ds)

		wₙ = wₙ ./ sum(abs.(wₙ),dims=2)
		wₙ = dropdims(mean(abs.(wₙ),dims=1),dims=1)
		prcp = dropdims(mean(prcp,dims=1),dims=1)

		ii = prcp .< (prcpRCEμ_T / 0.9)
		jj = prcp .> (prcpRCEμ_T * 0.9)

		for imode = 1 : nmode
			if !iszero(length(ii))
				anm[10].scatter(-mean(wₙ[imode,ii])./mean(wₙ[1,ii]),conDGW,c=blues_WTG[imode+1],s=10)
			end
			if !iszero(length(jj))
				anm[10].scatter(mean(wₙ[imode,jj])./mean(wₙ[1,jj]),conDGW,c=blues_WTG[imode+1],s=10)
			end
		end

	end

	anm[1].format(ylabel=L"$\tau$ / hr",ltitle="WTG")
	anm[2].format(ltitle="SWTG")
	anm[3].format(ltitle="DGW (Blossey)")
	anm[4].format(ltitle="DGW (Kuang)")
	anm[5].format(ylabel=L"$\alpha$",ltitle=L"DGW ($\partial_t\neq0$)")
	anm[6].format(ylabel=L"$\tau$ / hr")
	anm[10].format(ylabel=L"$\alpha$")

	for ax in anm
		ax.format(
			xlim=(-1.1,1.1),yscale="log",
			lrtitle="Wet",lltitle="Dry",urtitle="Wet",ultitle="Dry",
		)
	end

	for ii = vcat(1:2,6:7)
		anm[ii].format(ylim=(0.05,200))
	end
	for ii = vcat(3:5,8:10)
		anm[ii].format(ylim=(0.004,2500),ytickloc="r")
	end
	for ii = vcat(2:4,7:9)
		anm[ii].format(yticklabels=["","","",""],ytickminor=10:10:100,)
	end
	
	fnm.savefig(plotsdir("12-wn-norm.png"),transparent=false,dpi=400)
	load(plotsdir("12-wn-norm.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╠═c0205377-e7c3-4272-8dc0-68230da2118b
# ╟─4a85fc54-279c-449d-b84b-1dcb00d1ca74
# ╟─f032247e-7959-4a7a-8de2-7e09cf918b48
