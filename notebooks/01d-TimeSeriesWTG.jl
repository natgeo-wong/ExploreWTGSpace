### A Pluto.jl notebook ###
# v0.16.1

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
	using LinearAlgebra
	using NCDatasets
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the TroPrecLS project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 1d. Time-Series of 3D Variables

Here, we are doing an exploration of the time-series for some of the 3D variables, especially those that are important in solving for the WTG vertical velocity.
"

# ╔═╡ b6892634-9199-11eb-38d5-8dda8da16ed7
md"
### A. Just Doing a Double-Check

We recall the momentum damping equation (assuming $a_m$ is constant with height):

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Here, we just want to do a double check to see if the perturbation in $T_v$ can be directly inverted and solved to get $\omega'$, or if extra levels are needed.
"

# ╔═╡ d3b025e0-5b35-11eb-330a-5fbb2204da63
expname = "S1284km300"

# ╔═╡ 58095159-99b3-4810-b917-17be2093b8c2
config = "damping032"

# ╔═╡ 074370dd-64bf-4477-bde7-af1860e94097
imember = 1

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
begin
	dwtg = NCDataset(outstatname(expname,config,false,true,imember))
	time     = dwtg["time"][:] .- 80.5
	vert_WTG = dwtg["z"][:]
	plvl_WTG = dwtg["p"][:]
	pre_WTG  = dwtg["PRES"][:]
	wwtg_WTG = dwtg["WWTG"][:]
	tabs_WTG = dwtg["TABS"][:]
	qvap_WTG = dwtg["QV"][:]
	qcon_WTG = dwtg["QCOND"][:]
	close(dwtg)
end

# ╔═╡ 3ae60162-7009-43d8-af77-df987e2d7792
begin
	drce = NCDataset(outstatname("Control","S1284km300",false,false))
	vert_RCE = drce["z"][:]
	plvl_RCE = drce["p"][:]
	tabs_RCE = drce["TABS"][:]
	qvap_RCE = drce["QV"][:]
	qcon_RCE = drce["QCOND"][:]
	close(drce)
end

# ╔═╡ 0847d519-3f79-437b-97d8-bda02218fcd9
begin
	tvrt_WTG = tabs_WTG .* (1 .+ 0.61 * qvap_WTG/1000 - qcon_WTG/1000)
	tvrt_RCE = tabs_RCE .* (1 .+ 0.61 * qvap_RCE/1000 - qcon_RCE/1000)
	tvrt_RCE = dropdims(mean(tvrt_RCE[:,501:1000],dims=2),dims=2)
end

# ╔═╡ a2b66d2c-e038-4750-b4d2-28f20e6768ef
begin
	pplt.close(); fig,axs = pplt.subplots(nrows=3,aspect=4,axwidth=6)
	
	c1 = axs[1].contourf(
		time,plvl_WTG,wwtg_WTG,
		levels=vcat(-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50)/100,
		extend="both",cmap="RdBu_r"
	)
	axs[1].colorbar(c1,loc="r",label=L"$w_{wtg}$ / m s$^{-1}$")
	
	c1 = axs[2].contourf(
		time,plvl_WTG,tvrt_WTG .- tvrt_RCE,
		levels=vcat(-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50),
		extend="both",cmap="RdBu_r"
	)
	axs[2].colorbar(c1,loc="r",label=L"$T_v'$ / K")
	
	c1 = axs[3].contourf(
		time,plvl_WTG,pre_WTG .- plvl_WTG,
		levels=vcat(-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50),
		extend="both",cmap="RdBu_r"
	)
	axs[3].colorbar(c1,loc="r",label="p' / Pa")
	
	for ax in axs
		ax.format(yscale="log")
	end
	
	fig.savefig("test.png",transparent=false,dpi=150)
	load("test.png")
end

# ╔═╡ f5b70c65-431c-45cf-bd8e-c5ba078d76a1
md"### B. Try taking an inverse!"

# ╔═╡ 96b1d658-bb22-4e47-9102-8e749e508833
begin
	nlvl = length(plvl_WTG)
	Amat = Tridiagonal(ones(nlvl-3),ones(nlvl-2)*-2,ones(nlvl-3))
end

# ╔═╡ 138a0bf9-ae3f-45f9-98b8-ba0a7c9082f0
f = (tvrt_WTG .- tvrt_RCE) ./ pre_WTG

# ╔═╡ 41c033e2-ea1e-45ca-b646-edfa446c5adb
plvl_WTG

# ╔═╡ 35170d85-412a-499c-8e00-c63da6784d05
plvl_WTG .- plvl_RCE

# ╔═╡ b3a56224-8d6b-4456-ae0b-62be57381c89
u = zeros(nlvl,139)

# ╔═╡ 58be642c-789d-4283-bd52-872b86730c09
u[2:27,:] = Amat \ f[2:27,:]

# ╔═╡ 22eda698-d45b-4651-9270-3f0830465276
begin
	pplt.close(); f2,a2 = pplt.subplots(nrows=1,aspect=4,axwidth=6)
	
	c2 = a2[1].contourf(
		time,plvl_WTG,u,
		levels=vcat(-50,-20,-10,-5,-2,-1,0,1,2,5,10,20,50)/10,
		extend="both",cmap="RdBu_r"
	)
	a2[1].colorbar(c1,loc="r",label=L"Pa s$^{-1}$")
	# a2[1].format(yscale="log")
	
	f2.savefig("test2.png",transparent=false,dpi=150)
	load("test2.png")
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╠═d3b025e0-5b35-11eb-330a-5fbb2204da63
# ╠═58095159-99b3-4810-b917-17be2093b8c2
# ╠═074370dd-64bf-4477-bde7-af1860e94097
# ╠═55230f4a-7661-11eb-1c37-8b022b95e08e
# ╠═3ae60162-7009-43d8-af77-df987e2d7792
# ╠═0847d519-3f79-437b-97d8-bda02218fcd9
# ╠═a2b66d2c-e038-4750-b4d2-28f20e6768ef
# ╟─f5b70c65-431c-45cf-bd8e-c5ba078d76a1
# ╠═96b1d658-bb22-4e47-9102-8e749e508833
# ╠═138a0bf9-ae3f-45f9-98b8-ba0a7c9082f0
# ╠═41c033e2-ea1e-45ca-b646-edfa446c5adb
# ╠═35170d85-412a-499c-8e00-c63da6784d05
# ╠═b3a56224-8d6b-4456-ae0b-62be57381c89
# ╠═58be642c-789d-4283-bd52-872b86730c09
# ╠═22eda698-d45b-4651-9270-3f0830465276
