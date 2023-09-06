### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 77f51599-c5eb-4062-b621-3ba36f9549c3
begin
	using Pkg; Pkg.activate()
	using DrWatson
	
md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 984dd824-c22e-4ea4-bea2-4df5dc227333
begin
	@quickactivate "ExploreWTGSpace"
	using NCDatasets
	using NumericalIntegration
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ 35da8c94-3f5f-11ed-2ac8-e3fe7e417233
md"
# 04a. Constructing the CVR Curve
"

# ╔═╡ bbd4caf3-e8ad-4a7a-833d-889b2e1e5ade
begin
	configSST = vcat(299.86:0.01:299.9,300,301,302,305)
	nSST = length(configSST)
	md"Loading SST configurations ..."
end

# ╔═╡ 0ac4319c-b469-4aa3-b97e-eb06e362d096
function tair2qsat(T,P)

    tb = T - 273.15
    if tb <= 0
    	esat = exp(43.494 - 6545.8/(tb+278)) / (tb+868)^2
    else
    	esat = exp(34.494 - 4924.99/(tb+237.1)) / (tb+105)^1.57
    end


    r = 0.622 * esat / max(esat,P-esat)
    return r / (1+r)

end

# ╔═╡ 52926d36-1a13-4ea2-b675-5f19ebaee57b
begin
	pplt.close(); fig,axs = pplt.subplots(aspect=2,axwidth=3)
	
	for isst in configSST
	
		sstname = @sprintf("%6.2f",isst)
		sstname = replace(sstname,"."=>"d")
		sstname = "WTG_ExploreWTGSpace-SST$(sstname)K.nc"
	
		ds  = NCDataset(datadir("CVR","OUT_STAT",sstname))
		nt  = ds.dim["time"] - 1500
		swv = zeros(nt)
		pwv = ds["PW"][1501:end]
		prc = ds["PREC"][1501:end] / 24
		qv  = ds["QV"][:,1501:end] / 1000
		ta  = ds["TABS"][:,1501:end]
		pa  = ds["PRES"][:,1501:end] * 100
		qs  = tair2qsat.(ta,pa)
		mqs = mean(qs)
		for it = 1 : nt
			swv[it] = integrate(reverse(pa[:,it]),reverse(qs[:,it])) / 9.81 / 1000
		end
		swv = swv * 1000
		crh = mean(pwv ./ swv) * 100
		close(ds)
	
		axs[1].scatter(crh,mean(prc))
		
	end
	
	axs[1].format(ylim=(0.02,5),yscale="log",xlim=(70,100))
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ Cell order:
# ╟─35da8c94-3f5f-11ed-2ac8-e3fe7e417233
# ╟─77f51599-c5eb-4062-b621-3ba36f9549c3
# ╠═984dd824-c22e-4ea4-bea2-4df5dc227333
# ╟─bbd4caf3-e8ad-4a7a-833d-889b2e1e5ade
# ╠═0ac4319c-b469-4aa3-b97e-eb06e362d096
# ╟─52926d36-1a13-4ea2-b675-5f19ebaee57b
