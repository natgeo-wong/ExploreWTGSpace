### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 696cb854-ef36-11eb-2a3f-ef59f539d240
begin
	using Pkg; Pkg.activate()
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ c156225b-a907-47bd-8465-06e6fd3d8c09
begin
	@quickactivate "ExploreWTGSpace"
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ a5fc15fc-e8e3-47db-aa4b-32d73d57eddb
configlist = [
	"damping001","damping002","damping004","damping008","damping016",
	"damping032","damping064","damping128","damping256","damping512",
]

# ╔═╡ d975a956-a5f9-4096-8800-d7f10e9e955d
begin
	ncon = length(configlist)
	blu  = pplt.Colors("blues",(ncon+2))
	grn  = pplt.Colors("greens",(ncon+2))
	red  = pplt.Colors("reds",(ncon+2))
	ppl  = pplt.Colors("purples",(ncon+2))
	brn  = pplt.Colors("brown",(ncon+2))
	ylw  = pplt.Colors("fire_r",(ncon+2))
	orn  = pplt.Colors("orange2",(ncon+2))
	lgd  = Dict("frame"=>false,"ncols"=>1)
end

# ╔═╡ 6ec2a114-4362-46a3-aead-b858002d91f5
function plotcsf(expname,configlist,clrs,axsii)
	
	for ic in 1 : ncon
		imem = 0
		while imem < 100; imem += 1
			fnc = outstatname(expname,configlist[ic],false,true,imem)
			if isfile(fnc)
				_,p,t = retrievedims(fnc); t = t .- 80
				pr = retrievevar("PREC",fnc)./24
				pw = retrievevar("PW",fnc)
				ta = retrievevar("TABS",fnc)
				qv = retrievevar("QV",fnc)
				rh = calcrh(qv,ta,p)
				sw = calcswp(rh,qv,p)
				cr = pw ./ sw / 10
				if isone(imem) && (ic==(ncon-2))
					axsii.scatter(
						mean(cr[(end-99):end]),
						mean(pr[(end-99):end]),
						color=clrs[ic+1],s=10,
						label="$expname",legend="r",legend_kw=lgd
					)
				else
					axsii.scatter(
						mean(cr[(end-99):end]),
						mean(pr[(end-99):end]),
						color=clrs[ic+1],s=10
					)
				end
			end
		end
		
	end
	
end	

# ╔═╡ 26d3b6cb-193f-40d9-8b65-0185916f27cf
begin
	pplt.close()
	fts,ats = pplt.subplots(ncols=1,aspect=2,axwidth=4,sharex=0)
	
	plotcsf("P064km301d7",configlist,blu,ats[1])
	plotcsf("D064km301d7",configlist,grn,ats[1])
	plotcsf("H064km301d7",configlist,ylw,ats[1])
	plotcsf("T064km301d7",configlist,orn,ats[1])
	plotcsf("P064km305d0",configlist,red,ats[1])
	plotcsf("P064km295d0",configlist,ppl,ats[1])
	plotcsf("P128km301d7",configlist,brn,ats[1])
	ats[1].format(
		xlim=(30,100),xlabel="Column Relative Humidity / %",
		ylim=10. .^(-4.5,1.5),ylabel=L"Precipitation Rate / mm hr$^{-1}$",
		yscale="log"
	)
	
	fts.savefig(plotsdir("csfvprcp.png"),transparent=false,dpi=200)
	load(plotsdir("csfvprcp.png"))
end

# ╔═╡ Cell order:
# ╟─696cb854-ef36-11eb-2a3f-ef59f539d240
# ╟─c156225b-a907-47bd-8465-06e6fd3d8c09
# ╠═a5fc15fc-e8e3-47db-aa4b-32d73d57eddb
# ╠═d975a956-a5f9-4096-8800-d7f10e9e955d
# ╠═6ec2a114-4362-46a3-aead-b858002d91f5
# ╟─26d3b6cb-193f-40d9-8b65-0185916f27cf
