### A Pluto.jl notebook ###
# v0.14.2

using Markdown
using InteractiveUtils

# ╔═╡ 681658b0-5914-11eb-0d65-bbace277d145
begin
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 6dce35fc-5914-11eb-0ce2-0d4e164e1898
begin
	@quickactivate "ExploreWTGSpace"
	using Printf
	using Statistics
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 6b. Transitioning from RCE to WTG

In this notebook, we investigate and develop a way to implement the WTG forcing gradually in the System of Atmospheric Modelling.  The sudden introduction of WTG large-scale forcing often causes a model to enter a \"shocked\" state that unnaturally forces the model into a different state.  Here, we develop a method that gradually increases the strength of the WTG momentum-damping parameter from a near/pseudo-RCE state.
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
expname = "P128km301d7"

# ╔═╡ a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
begin
	configvec = [
		"damping001",
		"damping002",
		"damping004",
		"damping008",
		"damping016",
		"damping032",
		"damping064",
		"damping128",
		"damping256",
		"damping512",
	]
	ncon = length(configvec)
	blues = pplt.Colors("Blues",(ncon+2))
	reds  = pplt.Colors("Reds",(ncon+2))
	teals = pplt.Colors("Greens",(ncon+2))
	lgd = Dict("frame"=>false,"ncols"=>1)
md"Loading time dimension and defining the damping experiments ..."
end

# ╔═╡ 223b4286-8811-11eb-0e67-4da65e1999a5
function daymean(data)
	
	return dropdims(mean(reshape(data,1,:),dims=1),dims=1)
	
end

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	pplt.close()
	fts,ats = pplt.subplots(nrows=3,aspect=3,axwidth=4,hspace=0.2,sharey=0)
	
	for ic in 1 : ncon
		config = configvec[ic]
		config = replace(config,"damping"=>"")
		config = replace(config,"d"=>".")
		config = parse(Float64,config)
		imem = 0
		
		while imem < 100; imem += 1
			fnc = outstatname(expname,configvec[ic],false,true,imem)
			if isfile(fnc)
				_,_,t = retrievedims(fnc); t = t .- floor(t[1])
				td = daymean(t)
				pr = retrievevar("PREC",fnc)
				pa = retrievevar("AREAPREC",fnc)
				pa = daymean(pr ./ pa)
				pr = daymean(pr)
				sw = daymean(retrievevar("SWNS",fnc))
				lw = daymean(retrievevar("LWNS",fnc))
				sh = daymean(retrievevar("SHF",fnc))
				lh = daymean(retrievevar("LHF",fnc))
				pw = daymean(retrievevar("PW",fnc))
				seb = sw .- lw .- sh .- lh
				ats[1].plot(td,pr,lw=1,color=blues[ic+1])
				ats[3].plot(td,seb,lw=1,color=reds[ic+1])
				if imem == 1
					constr = @sprintf("%d",config)
					ats[2].plot(
						td,pw,lw=1,color=teals[ic+1],
						label=(L"$a_m =$" * " $(constr)"),
						legend="r",legend_kw=lgd
					)
				else
					ats[2].plot(td,pw,lw=1,color=teals[ic+1])
				end
			end
		end
		
	end
	
	ats[1].format(
		ylabel=L"Rainfall / mm day$^{-1}$",yscale="log",
		ylocator=10. .^(-3:3),
		# yscale_kw=Dict("linthresh"=>0.1),
		ylim=(0.005,200),
		suptitle=expname,ultitle="(a)"
	)
	
	ats[2].format(
		ylim=(0,75),ylabel="PW / mm",
		suptitle=expname,ultitle="(b)"
	)
	
	ats[3].format(
		xlim=(00,500),xlabel="Time / Days",
		ylim=(-250,250),ylabel=L"SEB / W m$^{-2}$",ylocator=(-3:3)*100,
		suptitle=expname,ultitle="(c)"
	)
	
	fts.savefig(plotsdir(
		"rce2wtg-$(expname).png"),
		transparent=false,dpi=150
	)
	load(plotsdir("rce2wtg-$(expname).png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╟─a63de98c-5b35-11eb-0a8f-b7a1ebd441b6
# ╟─223b4286-8811-11eb-0e67-4da65e1999a5
# ╠═55230f4a-7661-11eb-1c37-8b022b95e08e
