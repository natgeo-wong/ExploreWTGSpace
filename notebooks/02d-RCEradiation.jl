### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ c55e0539-a750-4205-8f36-a5ac8b1e0106
begin
	using Pkg; Pkg.activate()
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ cafa7b1e-1aff-4e90-bbaf-5076ae5d8416
begin
	@quickactivate "ExploreWTGSpace"
	using Statistics
	using PlutoUI
	using Printf
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	include(srcdir("samrad.jl"))
	
md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ 1722c69c-b534-11eb-09f3-4125b0813eab
md"
# 2d. Calculating Radiative Tendencies

This notebook will be used to calculate radiative tendencies from equilibrium profiles, to create a `rad` file to be input into the model such that non-interactive radiation can be ran.
"

# ╔═╡ 81f8007e-c17e-4b8e-b0a7-004c40ddd7d4
md"
### A. Definining the Model Configuration

There are two broad model configuration categories: (P)erpetual (INSOL)ation, and (D)iurnal (INSOL)ation.
"

# ╔═╡ f707eeb6-4948-4a09-9f1f-b8c118879ced
md"Is Diurnal? $(@bind isdiurnal PlutoUI.Slider(0:1))"

# ╔═╡ 446cfac8-26b7-46f9-a4f0-6bf0ec1ca710
md"Toggle Domain Size: $(@bind islarge PlutoUI.Slider(0:1))"

# ╔═╡ 2022819a-4f35-409d-9f44-499e0a174d06
md"Sea Surface Temperature: $(@bind issst PlutoUI.Slider(0:2))"

# ╔═╡ a8d202d7-d024-4217-a253-e2f42e5e38e8
begin
	if isone(isdiurnal)
		  insol = "D"
	else; insol = "P"
	end
	
	if isone(islarge)
		  domsize = "128"
	else; domsize = "064"
	end
	
	if iszero(issst)
		domsst = "295d0"
	elseif isone(issst)
		domsst = "301d7"
	else
		domsst = "305d0"
	end
		
	
	config = "$(insol)$(domsize)km$(domsst)"
end

# ╔═╡ f30ec381-6093-4cd0-a5c0-eac3bc616206
begin
	lvls = vcat(-10,-7.07,-5,-3.16,-2,-1.41,-1,-0.5,0.5,1,1.41,2,3.16,5,7.07,10)
	md"Defining universal plotting variables"
end

# ╔═╡ 6ff3200e-f014-46df-89e5-a8400eae0e1c
nen = 10

# ╔═╡ 186f97ce-20c3-4277-8284-c256386b4738
md"
### B. Extracting the Radiative Tendencies

From this initial spinup, I ran a $(nen)-member RCE ensemble, with initial conditions perturbed as per SAM's ensemble initialization (which allows for up to 100 members).  We use a model ensemble in order to account for the chaotic nature of convection and how that might perturb the final mean state.  Each member is run for 1000-1500 days, and we take the radiative tendencies of the last 250-500 days of the model run, and create `rad` forcing file that overrides the old file.
"

# ╔═╡ 1551de54-c09b-44f3-bcd1-15d8189801a8
nendays = 500

# ╔═╡ 5594019f-c02a-48ee-b0f0-f39a6515901e
begin
	z_en,_,t_en = retrievedims("Control","$(config)",isensemble=true,member=1)
	nz_en = length(z_en); nt_en = length(t_en)
	pre_en = zeros(nz_en,nt_en,nen); plevel = zeros(nz_en,nen)
	rad_en = zeros(nz_en,nt_en,nen)
	for imem = 1 : nen
		rad_en[:,:,imem] = retrievevar("RADQR","Control","$(config)",isensemble=true,member=imem)
		pre_en[:,:,imem] = retrievevar("PRES","Control","$(config)",isensemble=true,member=imem)
		plevel[:,imem] = retrievevar("p","Control","$(config)",isensemble=true,member=imem)[:,end]
	end
	
	plevel = dropdims(mean(plevel,dims=2),dims=2)
md"Loading data from the $(nen)-member ensemble for the $(config) RCE run ..."
end

# ╔═╡ a8560cec-743b-4138-88ed-501dbbd5ad3a
begin
	rtend_en = dropdims(mean(rad_en[:,(end-nendays+1):end,:],dims=2),dims=2)
	radts_en = dropdims(mean(rad_en,dims=3),dims=3)
	pts_en   = dropdims(mean(pre_en[:,(end-nendays+1):end,:],dims=2),dims=2)
md"Calculating the radiative tendencies from the last $nendays days of model run ..."
end

# ╔═╡ 193e20a8-7469-4583-915c-a4ccbebf82d8
begin
	
	pplt.close(); fen,aen = pplt.subplots(aspect=2.25,axwidth=4,sharex=0)
	
	crad = aen[1].contourf(
		t_en.-80,plevel,radts_en,
		# t.-80,p,tem .- mean(tem[:,(end-100+1):end],dims=2),
		# t[5:(end-4)].-80,p,temn .- tob[:,1],
		cmap="RdBu_r",cmap_kw=Dict("alpha"=>(1,1,1,1,1,1,0,1,1,1,1,1,1)),
		extend="both",
		levels=lvls/10
	)
	aen[1].format(
		xlim=(0,nt_en),
		ylim=(1010,25),yscale="log",
		ylabel="Pressure / hPa",
		ultitle="(b)"
	)
	
	pen = aen[1].panel("l",space=0.1,width="5em")
	
	for im = 1 : nen
		pen.scatter(rtend_en[:,im],pts_en[:,im],s=2,c="gray")
	end
	
	pen.plot(
		dropdims(mean(rtend_en,dims=2),dims=2),
		dropdims(mean(pts_en,dims=2),dims=2),c="k"
	)
	pen.format(
		xlim=(-1.5,1.5),xlocator=(-2:2),xlabel=L"T - T$_{OBS}$ / K",
		ylim=(1010,25),yscale="log",ylabel="Pressure / hPa",
		suptitle="Model Ensemble Equilibrium RCE | $config",ultitle="(a)"
	)
	
	aen[1].colorbar(crad,loc="r",width=0.2,locator=[-10,-5,-2,-1,0,1,2,5,10]/10)
	fen.savefig(plotsdir("rcetdts-$(config)-ensemble.png"),transparent=false,dpi=200)
	load(plotsdir("rcetdts-$(config)-ensemble.png"))
	
end

# ╔═╡ aade6085-6ee3-468b-93c0-a16265ec379e
begin
	rad_μ = dropdims(mean(rad_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	pre_μ = dropdims(mean(pre_en[:,(end-nendays+1):end,:],dims=(2,3)),dims=(2,3))
	
	raddata = zeros(nz_en,2); raddata[:,1] .= pre_μ; raddata[:,2] .= rad_μ/86400
end

# ╔═╡ f721d659-d9d7-4d5c-98b6-5f29e9a6d0cc
md"Create RAD file from ensemble? $(@bind dorad PlutoUI.Slider(0:1))"

# ╔═╡ 6ef336d3-fec3-40da-9ac1-07a7b66ff48a
if isone(dorad)
	  createradmean("$(config)",raddata)
	  md"Creating the RAD file $(config) from ensemble simulations ..."
else; md"We have decided not to create the ensemble RAD file $(config) yet ..."
end

# ╔═╡ Cell order:
# ╟─1722c69c-b534-11eb-09f3-4125b0813eab
# ╟─c55e0539-a750-4205-8f36-a5ac8b1e0106
# ╟─cafa7b1e-1aff-4e90-bbaf-5076ae5d8416
# ╟─81f8007e-c17e-4b8e-b0a7-004c40ddd7d4
# ╟─f707eeb6-4948-4a09-9f1f-b8c118879ced
# ╟─446cfac8-26b7-46f9-a4f0-6bf0ec1ca710
# ╟─2022819a-4f35-409d-9f44-499e0a174d06
# ╟─a8d202d7-d024-4217-a253-e2f42e5e38e8
# ╟─f30ec381-6093-4cd0-a5c0-eac3bc616206
# ╟─186f97ce-20c3-4277-8284-c256386b4738
# ╟─6ff3200e-f014-46df-89e5-a8400eae0e1c
# ╟─1551de54-c09b-44f3-bcd1-15d8189801a8
# ╟─5594019f-c02a-48ee-b0f0-f39a6515901e
# ╠═a8560cec-743b-4138-88ed-501dbbd5ad3a
# ╟─193e20a8-7469-4583-915c-a4ccbebf82d8
# ╟─aade6085-6ee3-468b-93c0-a16265ec379e
# ╟─f721d659-d9d7-4d5c-98b6-5f29e9a6d0cc
# ╟─6ef336d3-fec3-40da-9ac1-07a7b66ff48a
