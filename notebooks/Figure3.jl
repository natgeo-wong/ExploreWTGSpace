### A Pluto.jl notebook ###
# v0.19.27

using Markdown
using InteractiveUtils

# ╔═╡ 0ed8f41a-74fe-11ed-1f19-85b3e328530b
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson in order to ensure reproducibility between different machines ..."
end

# ╔═╡ 000a8d8d-24c7-40f6-a47a-8a299bb508f5
begin
	@quickactivate "ExploreWTGSpace"
	using DelimitedFiles
	using NCDatasets
	using StatsBase
	using Printf
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow

	pplt = pyimport("proplot")

	md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ 0dcabe70-3b22-4b76-a42e-80ab872da27d
md"
# Figure 3. Idealized Static Energies
"

# ╔═╡ 9b813dbb-80dc-4e52-9588-5b8ea360e483
zt = 0 : 0.01 : 1

# ╔═╡ 230c7c11-7f41-4bfa-9094-5d187b5cf741
zz = 0 : 0.01 : 1.1

# ╔═╡ d1cf9344-ef1f-4915-a674-f9e02b0cb7e0
id_mse = 0.5 .+(zz.^2 .-zz)./2 .- 0.375

# ╔═╡ 81106d9d-7141-4aeb-8539-2e704e031ed6
begin
	id_dse = deepcopy(id_mse)
	id_dse[zz.<0.5] .= 0 .- id_dse[zz.<0.5]
end

# ╔═╡ 97299df2-efb6-44c8-8603-c9d2b4de6f0a
id_gmse = (zz.-0.5)./10

# ╔═╡ a4e8272c-4a4d-4f4e-8ef9-2ce607a499fd
begin
	id_gdse = collect(deepcopy(id_gmse))
	id_gdse[zz.<0.5] .= -1 * id_gdse[zz.<0.5]
end

# ╔═╡ 00f4548d-c8cd-455e-92f0-9d3c92499211
begin
	w1_wet = 0.1*sin.(zz*pi);   w1_wet[zz.>1] .= 0
	w2_wet = 0.1*sin.(zz*2*pi); w2_wet[zz.>1] .= 0
	w1_dry = -0.07*sin.(zz*pi);   w1_dry[zz.>1] .= 0
	w2_dry = -0.07*sin.(zz*2*pi); w2_dry[zz.>1] .= 0
end

# ╔═╡ 6e71c910-87ba-49cf-bfc2-33d2c5323051
lgd = Dict("ncols"=>1,"frame"=>false)

# ╔═╡ 0dc0c36f-9bc3-4181-a552-dd9e097eb00d
begin
	pplt.close
	fig,axs = pplt.subplots(
		ncols=5,aspect=0.5,axwidth=1,
		sharey=0,sharex=0,wspace=[3,1,3,1]
	)
	
	axs[1].plot(c="blue6", alpha=0.5,w1_wet,zz,label=L"$W_1$ (Moist)",legend="b")
	axs[1].plot(c="yellow6",alpha=0.5,w2_dry,zz,label=L"$W_2$ (Stratiform)",legend="b",legend_kw=lgd)
	axs[1].format(xlocator=[0],ultitle="(a)")

	axs[2].plot(id_mse,zz,c="k",label="h",legend="b",legend_kw=lgd)
	axs[2].plot(id_dse,zz,c="gray",linestyle="--",label="s",legend="b")
	axs[2].format(xlocator=[],ultitle="(b1)")

	axs[3].plot(id_gmse,zz,c="k",label=L"\partial_zh",legend="b",legend_kw=lgd)
	axs[3].plot(id_gdse,zz,c="gray",linestyle="--",label=L"\partial_zs",legend="b")
	axs[3].format(xlocator=[0],ultitle="(b2)")

	axs[4].plot(id_gmse.*w1_wet.*25,zz,c="blue6",label=L"W_1 \cdot \partial_zh",legend="b",legend_kw=lgd)
	axs[4].plot(id_gdse.*w1_wet.*25,zz,c="blue2",linestyle="--",label=L"W_1 \cdot \partial_zs",legend="b")
	axs[4].format(xlocator=[0],ultitle="(c1)")

	axs[5].plot(id_gmse.*w2_dry.*25,zz,c="yellow6",label=L"W_2 \cdot \partial_zh",legend="b",legend_kw=lgd)
	axs[5].plot(id_gdse.*w2_dry.*25,zz,c="yellow2",linestyle="--",label=L"W_2 \cdot \partial_zs",legend="b")
	axs[5].format(xlocator=[0],ultitle="(c2)")

	for ax in axs
		ax.format(
			xlim=(-0.2,0.2),ylim=(0,1.2),xloc="bottom",grid=false,yloc="zero",
			ylocator=0:1,yminorlocator=[],yticklabels=["",L"z_t"]
		)
	end
	
	fig.savefig(projectdir("figures","fig3-idealisedstaticenergies.png"),transparent=false,dpi=400)
	load(projectdir("figures","fig3-idealisedstaticenergies.png"))
end

# ╔═╡ Cell order:
# ╟─0dcabe70-3b22-4b76-a42e-80ab872da27d
# ╟─0ed8f41a-74fe-11ed-1f19-85b3e328530b
# ╟─000a8d8d-24c7-40f6-a47a-8a299bb508f5
# ╠═9b813dbb-80dc-4e52-9588-5b8ea360e483
# ╠═230c7c11-7f41-4bfa-9094-5d187b5cf741
# ╠═d1cf9344-ef1f-4915-a674-f9e02b0cb7e0
# ╠═81106d9d-7141-4aeb-8539-2e704e031ed6
# ╠═97299df2-efb6-44c8-8603-c9d2b4de6f0a
# ╠═a4e8272c-4a4d-4f4e-8ef9-2ce607a499fd
# ╠═00f4548d-c8cd-455e-92f0-9d3c92499211
# ╠═6e71c910-87ba-49cf-bfc2-33d2c5323051
# ╟─0dc0c36f-9bc3-4181-a552-dd9e097eb00d
