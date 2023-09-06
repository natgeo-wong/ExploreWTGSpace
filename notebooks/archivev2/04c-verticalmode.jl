### A Pluto.jl notebook ###
# v0.19.16

using Markdown
using InteractiveUtils

# ╔═╡ 0ed8f41a-74fe-11ed-1f19-85b3e328530b
begin
	using Pkg; Pkg.activate()
	using DrWatson
end

# ╔═╡ 000a8d8d-24c7-40f6-a47a-8a299bb508f5
begin
	@quickactivate "AGU2022"
	using DelimitedFiles
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow

	pplt = pyimport("proplot")
end

# ╔═╡ ffc435cd-667e-44a9-b3eb-84443128aca5
begin
	z = 0 : 0.01 : 1
end

# ╔═╡ cc65882f-8725-4f20-a5fa-2d4fb3f75d7a
begin
	pplt.close
	f1,a1 = pplt.subplots([2,1,3],aspect=0.5,axwidth=0.75,sharey=0,wspace=5,sharex=0)

	a1[1].plot(c="k",alpha=0.3,sin.(z*pi),z)
	a1[2].plot(c="k",alpha=0.3,sin.(z*pi),z)
	a1[3].plot(c="k",alpha=0.3,sin.(z*pi),z,legend="r",label="1st Baroclinic Mode",legend_kw=Dict("ncols"=>1,"frame"=>false))

	a1[1].plot(c="blue6",alpha=0.4,sin.(z*2*pi),z)
	a1[2].plot(c="blue6",alpha=0.4,0.25*sin.(z*2*pi),z)
	a1[3].plot(c="blue6",alpha=0.4,sin.(z*2*pi),z,legend="r",label="2nd Baroclinic Mode")
	
	a1[1].plot(c="blue6",sin.(z*pi).+sin.(z*2*pi),z)
	a1[2].plot(c="blue6",sin.(z*pi).+sin.(z*2*pi)*0.25,z)
	a1[3].plot(c="blue6",sin.(z*pi).+sin.(z*2*pi),z,label="Total",legend="r")

	a1[1].format(xlabel="T'")
	a1[2].format(xlabel=L"\omega'")
	a1[3].format(xlabel=L"w'")

	for ax in a1
		ax.format(
			xlim=(-1.9,1.9),ylim=(0,1.2),yloc="zero",xloc="bottom",
			xlocator=-2:2,ylocator=0:1,grid=false,
			yminorlocator=[],yticklabels=["",L"z_t"]
		)
	end
	
	f1.savefig(plotsdir("verticalmode_congestus.png"),transparent=false,dpi=400)
	load(plotsdir("verticalmode_congestus.png"))
end

# ╔═╡ 33ae0e01-dc00-45db-93f5-6c7be26697cd
begin
	pplt.close
	f2,a2 = pplt.subplots([2,1,3],aspect=0.5,axwidth=0.75,sharey=0,wspace=5,sharex=0)

	a2[1].plot(c="k",alpha=0.3,sin.(z*pi),z)
	a2[2].plot(c="k",alpha=0.3,sin.(z*pi),z)
	a2[3].plot(c="k",alpha=0.3,sin.(z*pi),z,legend="r",label="1st Baroclinic Mode",legend_kw=Dict("ncols"=>1,"frame"=>false))

	a2[1].plot(c="yellow6",alpha=0.4,-sin.(z*2*pi),z)
	a2[2].plot(c="yellow6",alpha=0.4,-sin.(z*2*pi)*0.25,z)
	a2[3].plot(c="yellow6",alpha=0.4,-sin.(z*2*pi),z,legend="r",label="2nd Baroclinic Mode")
	
	a2[1].plot(c="yellow6",sin.(z*pi).-sin.(z*2*pi),z)
	a2[2].plot(c="yellow6",sin.(z*pi).-sin.(z*2*pi)*0.25,z)
	a2[3].plot(c="yellow6",sin.(z*pi).-sin.(z*2*pi),z,legend="r",label="Total")

	a2[1].format(xlabel="T'")
	a2[2].format(xlabel=L"\omega'")
	a2[3].format(xlabel=L"w'")

	for ax in a2
		ax.format(
			xlim=(-1.9,1.9),ylim=(0,1.2),yloc="zero",xloc="bottom",
			xlocator=-2:2,ylocator=0:1,grid=false,
			yminorlocator=[],yticklabels=["",L"z_t"]
		)
	end
	
	f2.savefig(plotsdir("verticalmode_stratiform.png"),transparent=false,dpi=400)
	load(plotsdir("verticalmode_stratiform.png"))
end

# ╔═╡ 9903c976-e59a-44b4-959d-37fee3157109


# ╔═╡ 0dc0c36f-9bc3-4181-a552-dd9e097eb00d


# ╔═╡ 7f8b325e-6f5b-432a-997d-7b040eda4692


# ╔═╡ Cell order:
# ╟─0ed8f41a-74fe-11ed-1f19-85b3e328530b
# ╠═000a8d8d-24c7-40f6-a47a-8a299bb508f5
# ╠═ffc435cd-667e-44a9-b3eb-84443128aca5
# ╠═cc65882f-8725-4f20-a5fa-2d4fb3f75d7a
# ╠═33ae0e01-dc00-45db-93f5-6c7be26697cd
# ╠═9903c976-e59a-44b4-959d-37fee3157109
# ╠═0dc0c36f-9bc3-4181-a552-dd9e097eb00d
# ╠═7f8b325e-6f5b-432a-997d-7b040eda4692
