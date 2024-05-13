### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ df810659-6173-483f-912c-523b022d641e
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 46faa412-5ade-11eb-3c37-23a7e59037a0
begin
	@quickactivate "ExploreWTGSpace"
	using Statistics
	using PlutoUI
	using Printf
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("samsnd.jl"))
	
md"Loading modules for the ExploreWTGSpace paper ..."
end

# ╔═╡ 9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
md"
# 1. Finding an Equilibrium RCE State

This notebook will be used to help find an equilibrium RCE state that roughly returns the initial sounding profile used to initialize the model.

SAM reads in \"snd\" (sounding profiles of humidity and potential temperature) files to initialize the model upon spinup.  We run models in RCE mode over 500-2000 day periods and obtain the profile of absolute temperature for the last `ndy` days.  We then compare it to the `TABSOBS` that was fed into the model, and find the temperature difference.  If the difference is too large, then we extract the vertical sounding profile, override the old \"snd\" file and rerun the model.
"

# ╔═╡ b2670c08-81e5-11eb-324e-2b923b289a04
md"
### A. Definining the Model Configuration

There are two broad model configuration categories: (P)erpetual (INSOL)ation, and (D)iurnal (INSOL)ation.
"

# ╔═╡ 505aec83-2f74-4f12-be40-b360cfb1d2a8
@bind prefix Select([
	"P" => "Perpetual Insolation (P)",
	"D" => "Diurnal Insolation (D)",
	"T" => "Temperature Tendency (T)",
])

# ╔═╡ 427511d0-88cb-11eb-2a40-019c91ee1401
if prefix == "D"
	md"Toggle Domain Size: $(@bind dsize PlutoUI.Slider(0:1,default=0))"
else
	md"Toggle Domain Size: $(@bind dsize PlutoUI.Slider(-1:0,default=0))"
end

# ╔═╡ 20d88eb4-8035-43d3-bfd1-084e38726b4b
if prefix == "P" && iszero(dsize)
	md"Toggle Horizontal Resolution: $(@bind hres PlutoUI.Slider(1:2,default=1))"
else
	hres = 1
	md""
end

# ╔═╡ 8b3016d6-bed8-4161-aa24-90db24b34920
if prefix == "P" && !isone(hres)
	md"Toggle Vertical Levels: $(@bind isvert PlutoUI.Slider(1:5,default=4))"
else
	isvert=0
	md""
end

# ╔═╡ 9cf41788-9587-4cdc-9761-93a9735ee10d
if (prefix == "P" || prefix == "T") && isone(hres) && iszero(dsize)
	md"Is 2D? $(@bind is2D PlutoUI.Slider(0:1,default=0))"
else
	is2D = 0
	md""
end

# ╔═╡ 6497001a-6400-4529-88d7-1a75597e5df5
if (prefix == "P" || prefix == "T") && isone(hres) && iszero(dsize) && iszero(is2D)
	md"Is Fixed Surface Flux? $(@bind isfsf PlutoUI.Slider(0:1,default=0))"
else
	isfsf = 0
	md""
end

# ╔═╡ 5940bed1-cf63-4caa-a954-b4aa8d3ab459
md"Domain Size = $(128*2. ^dsize) x $(128*2. ^dsize) km"

# ╔═╡ 8697387f-a384-41a6-bba1-f61f47056bd0
md"Horizontal Resolution: $(2. ^hres) km"

# ╔═╡ fab65122-bf01-44d7-8a92-e0b842983a98
begin
	if prefix == "P" && !isone(hres)
		vertvec = [28,42,54,64,144]
		nvert = vertvec[isvert]
	else
		nvert = 64
	end
	md"Number of Vertical Levels: $(nvert)"
end

# ╔═╡ ac8b9d4c-5ade-11eb-06f4-33bff063bbde
begin
	config = "$(prefix)$(@sprintf("%03d",128*2. ^dsize))$(@sprintf("%d",2. ^hres))km300V$(nvert)"
	if isone(is2D); config = "$(config)_2D" end
	if isone(isfsf); config = "$(config)_FSF" end
	nen = 10
	
md"**Experiment Set:** $config | **Number of Control Members**: $nen"
end

# ╔═╡ 7842e150-822b-11eb-1ded-f35ee4cc6d8c
begin
	arr = [[5,1,2,2,2,2,2,2],[6,3,4,4,4,4,4,4]]
	lvls = vcat(-5,-2,-1,-0.5,-0.2,-0.1,0.1,0.2,0.5,1,2,5)
	md"Defining universal plotting variables"
end

# ╔═╡ 99dca670-81e5-11eb-24b8-cf1b23d8c1f7
md"
### B. Running a $(nen)-member RCE ensemble

From this initial spinup, I ran a $(nen)-member RCE ensemble, with initial conditions perturbed as per SAM's ensemble initialization (which allows for up to 100 members).  We use a model ensemble in order to account for the chaotic nature of convection and how that might perturb the final mean state.  Each member is run for 1000-1500 days, and we take the temperature and humidity profiles of the last 250-500 days of the model run, and create a new sounding profile that overrides the old file.

Should the difference in the final averaged sounding profile be vastly different from the initial sounding profile, the model ensemble is run again, with the final averaged sounding profile as the new initial sounding profile.
"

# ╔═╡ 401a6a3a-8444-11eb-3ee8-594591ed5aa9
nendays = 1900

# ╔═╡ c658d650-8614-11eb-252c-c3066fb1d506
ptrop = 70

# ╔═╡ 7d905176-81e8-11eb-20d3-b9be287472f5
begin
	z_en,_,t_en = retrievedims("RCE","$(config)",isensemble=true,member=1)
	nz_en = length(z_en); nt_en = length(t_en);
	tbi_en = zeros(nz_en,nt_en,nen); tob_en = zeros(nz_en,nen)
	qbi_en = zeros(nz_en,nt_en,nen); qob_en = zeros(nz_en,nen)
	tem_en = zeros(nz_en,nt_en,nen)
	qvp_en = zeros(nz_en,nt_en,nen); rh_en  = zeros(nz_en,nt_en,nen)
	pre_en = zeros(nz_en,nt_en,nen); plevel = zeros(nz_en,nen)
	prc_en = zeros(nt_en,nen);
	shf_en = zeros(nt_en,nen);
	lhf_en = zeros(nt_en,nen);
	for imem = 1 : nen
		tbi_en[:,:,imem] = retrieve3Dvar("TBIAS","RCE","$(config)",isensemble=true,member=imem)[:,1:nt_en]
		qbi_en[:,:,imem] = retrieve3Dvar("QBIAS","RCE","$(config)",isensemble=true,member=imem)[:,1:nt_en]
		tem_en[:,:,imem] = retrieve3Dvar("TABS","RCE","$(config)",isensemble=true,member=imem)[:,1:nt_en]
		qvp_en[:,:,imem] = retrieve3Dvar("QV","RCE","$(config)",isensemble=true,member=imem)[:,1:nt_en]
		pre_en[:,:,imem] = retrieve3Dvar("PRES","RCE","$(config)",isensemble=true,member=imem)[:,1:nt_en]
		tob_en[:,imem] = retrieve3Dvar("TABSOBS","RCE","$(config)",isensemble=true,member=imem)[:,end]
		qob_en[:,imem] = retrieve3Dvar("QVOBS","RCE","$(config)",isensemble=true,member=imem)[:,end]
		plevel[:,imem] = retrieve2Dvar("p","RCE","$(config)",isensemble=true,member=imem)
		prc_en[:,imem] = retrieve2Dvar("PREC","RCE","$(config)",isensemble=true,member=imem)[1:nt_en]
		rh_en[:,:,imem]= calcrh(qvp_en[:,:,imem],tem_en[:,:,imem],plevel[:,imem])/10
		shf_en[:,imem] = retrieve2Dvar("SHF","RCE","$(config)",isensemble=true,member=imem)[1:nt_en]
		lhf_en[:,imem] = retrieve2Dvar("LHF","RCE","$(config)",isensemble=true,member=imem)[1:nt_en]
	end
	
	tob_en = dropdims(mean(tob_en,dims=2),dims=2)
	qob_en = dropdims(mean(qob_en,dims=2),dims=2)
	plevel = dropdims(mean(plevel,dims=2),dims=2)
md"Loading data from the $(nen)-member ensemble for the $(config) RCE run ..."
end

# ╔═╡ ad523b4e-81ee-11eb-2d10-a984d5983471
begin
	tdiff_en = dropdims(mean(tbi_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	tdts_en  = dropdims(mean(tbi_en,dims=3),dims=3)
	qdiff_en = dropdims(mean(qbi_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	qdiff_en = qdiff_en ./ qob_en*100
	qdts_en  = dropdims(mean(qbi_en ./ qob_en*100,dims=3),dims=3)
	pts_en   = dropdims(mean(pre_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	tabs_en  = dropdims(mean(tem_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	relh_en  = dropdims(mean(rh_en[:,(end-nendays+1):end,:],dims=2),dims=2)
md"Calculating difference between `TABS` and `TABSOBS` ..."
end

# ╔═╡ f35ab14c-8226-11eb-29b4-c3b7130f3733
begin
	ip_en = plevel .> ptrop
	trms_en = sqrt(mean(tdiff_en[ip_en,:].^2)); trms_en = @sprintf("%.3f",trms_en)
	qrms_en = sqrt(mean(qdiff_en[ip_en,:].^2)); qrms_en = @sprintf("%.3f",qrms_en)
	
md"Assuming that the tropopause in the tropics can sometimes reach as high as $ptrop hPa, the root-mean-square of the temperature difference between the model temperature `TABS` and the observed temperature `TABSOBS` below the tropopause is $(trms_en) K (i.e., we take calculate the root-mean-square using only levels where p > $(ptrop) hPa).  The profile of the temperature difference is shown below:"
end

# ╔═╡ a3650e8a-822b-11eb-29ca-cd01b90e8099
begin
	
	pplt.close()
	fen,aen = pplt.subplots(arr,aspect=1/3,axwidth=0.6,sharex=0,wspace=1)
	
	for im = 1 : nen
		aen[1].scatter(tdiff_en[:,im],pts_en[:,im],s=1,c="gray")
	end
	
	aen[1].plot(
		dropdims(mean(tdiff_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[1].format(
		xlim=(-1.5,1.5)./10,xlocator=(-2:2)./10,xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",
		suptitle="Model Ensemble Equilibrium RCE | $config",ultitle="(a)"
	)
	
	cten = aen[2].pcolormesh(
		t_en.-t_en[1],plevel,tdts_en,
		# t.-80,p,tem .- mean(tem[:,(end-100+1):end],dims=2),
		# t[5:(end-4)].-80,p,temn .- tob[:,1],
		cmap="RdBu_r",
		extend="both",
		levels=lvls/10
	)
	aen[2].format(
		xlim=(0,nt_en),
		ylim=(1010,25),yscale="log",
		ylabel="Pressure / hPa",
		urtitle=L"T$_{RMS}$" * " = $(trms_en) K",ultitle="(b)"
	)
	
	for im = 1 : nen
		aen[3].scatter(qdiff_en[:,im],pts_en[:,im],s=1,c="gray")
	end
	
	aen[3].plot(
		dropdims(mean(qdiff_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[3].format(
		xlim=(-1.5,1.5),xlocator=(-2:2),
		xlabel=L"qr = $\frac{q - q_{OBS}}{q_{OBS}}$ / %",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",ultitle="(c)"
	)
	
	cqen = aen[4].pcolormesh(
		t_en.-t_en[1],plevel,qdts_en,
		# t.-80,p,qvp .- mean(qvp[:,(end-100+1):end],dims=2),
		cmap="drywet",
		extend="both",
		levels=lvls
	)
	aen[4].format(
		xlim=(0,nt_en),
		ylim=(1010,25),yscale="log",
		xlabel="time / days",ylabel="Pressure / hPa",
		urtitle=L"$qr_{RMS}$" * " = $(qrms_en) %",ultitle="(d)"
	)
	
	for im = 1 : nen
		aen[5].scatter(tabs_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	aen[5].plot(
		dropdims(mean(tabs_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[5].format(xlim=(180,320),xlabel="T / K")
	
	for im = 1 : nen
		aen[6].scatter(relh_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	aen[6].plot(
		dropdims(mean(relh_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	aen[6].format(xlim=(0,120),xlabel="r / %")
	
	aen[2].colorbar(cten,loc="r",width=0.2)
	aen[4].colorbar(cqen,loc="r",width=0.2)
	fen.savefig(plotsdir("01a-ensemble-$(config).png"),transparent=false,dpi=400)
	load(plotsdir("01a-ensemble-$(config).png"))
	
end

# ╔═╡ 11913f26-8e7c-11eb-1b37-0f9c8e7b7331
md"The domain mean precipitation is $(mean(prc_en)) mm/day"

# ╔═╡ caf84abe-581d-4859-baf5-81c75d133fd4
begin
	shfstring = @sprintf("%.3f",mean(shf_en))
	lhfstring = @sprintf("%.3f",mean(lhf_en))
	md"The domain mean Surface and Latent Heat Fluxes are $(shfstring) W m^-2 and $(lhfstring) W m^-2 respectively"
end

# ╔═╡ 4db59bf0-82ec-11eb-0374-81a982a74216
md"
### D. Creating Sounding from Ensemble
"

# ╔═╡ 549c7744-82ed-11eb-03e8-7bb72c725856
begin
	qvp_μ = dropdims(mean(qvp_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	tem_μ = dropdims(mean(tem_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	pre_μ = dropdims(mean(pre_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	
	pot_μ = tem_μ .* (1000 ./pre_μ).^(287/1004)
	
	snddata = zeros(nz_en,6)
	snddata[:,1] .= z_en;  snddata[:,2] .= pre_μ
	snddata[:,3] .= pot_μ; snddata[:,4] .= qvp_μ
	md"Creating sounding data from ensemble data ..."
end

# ╔═╡ eba0fc7a-82eb-11eb-104e-d17de6c6c0de
md"Create SND file from ensemble? $(@bind dosnden PlutoUI.Slider(0:1))"

# ╔═╡ f53da85a-82eb-11eb-0f42-ad74458f849f
if isone(dosnden)
	  createsndmean("$(config)",snddata,psfc=1009.32)
	  md"Creating the SND file $(config) from ensemble simulations ..."
else; md"We have decided not to create the ensemble SND file $(config) yet ..."
end

# ╔═╡ Cell order:
# ╟─9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
# ╟─df810659-6173-483f-912c-523b022d641e
# ╟─46faa412-5ade-11eb-3c37-23a7e59037a0
# ╟─b2670c08-81e5-11eb-324e-2b923b289a04
# ╟─505aec83-2f74-4f12-be40-b360cfb1d2a8
# ╟─427511d0-88cb-11eb-2a40-019c91ee1401
# ╟─20d88eb4-8035-43d3-bfd1-084e38726b4b
# ╟─8b3016d6-bed8-4161-aa24-90db24b34920
# ╟─9cf41788-9587-4cdc-9761-93a9735ee10d
# ╟─6497001a-6400-4529-88d7-1a75597e5df5
# ╟─5940bed1-cf63-4caa-a954-b4aa8d3ab459
# ╟─8697387f-a384-41a6-bba1-f61f47056bd0
# ╟─fab65122-bf01-44d7-8a92-e0b842983a98
# ╟─ac8b9d4c-5ade-11eb-06f4-33bff063bbde
# ╟─7842e150-822b-11eb-1ded-f35ee4cc6d8c
# ╟─99dca670-81e5-11eb-24b8-cf1b23d8c1f7
# ╠═401a6a3a-8444-11eb-3ee8-594591ed5aa9
# ╠═c658d650-8614-11eb-252c-c3066fb1d506
# ╟─7d905176-81e8-11eb-20d3-b9be287472f5
# ╟─ad523b4e-81ee-11eb-2d10-a984d5983471
# ╟─f35ab14c-8226-11eb-29b4-c3b7130f3733
# ╟─a3650e8a-822b-11eb-29ca-cd01b90e8099
# ╟─11913f26-8e7c-11eb-1b37-0f9c8e7b7331
# ╟─caf84abe-581d-4859-baf5-81c75d133fd4
# ╟─4db59bf0-82ec-11eb-0374-81a982a74216
# ╟─549c7744-82ed-11eb-03e8-7bb72c725856
# ╟─eba0fc7a-82eb-11eb-104e-d17de6c6c0de
# ╟─f53da85a-82eb-11eb-0f42-ad74458f849f
