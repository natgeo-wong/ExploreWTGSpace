### A Pluto.jl notebook ###
# v0.15.1

using Markdown
using InteractiveUtils

# ╔═╡ 3f630f34-1d64-11ec-073c-316ff2a18be2
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to maintain reproducibility of project ..."
end

# ╔═╡ 5e772289-f828-4876-8456-edc9824716c6
begin
	@quickactivate "ExploreWTGSpace"
	using NCDatasets
	using Statistics
	
	using PyCall, LaTeXStrings
	using PNGFiles, ImageShow
	pplt = pyimport("proplot")
	
	md"Loading modules for the SelfAggregation Project ..."
end

# ╔═╡ 5647a8b2-ab77-4cab-b1e0-3d94f0c6d332
begin
	fol0Dl = datadir("RCEMIP","large","0D")
	fol1Dl = datadir("RCEMIP","large","1D")
	fol2Dl = datadir("RCEMIP","large","2D")
	fol3Dl = datadir("RCEMIP","large","3D")
	md"Folders containing RCEMIP data for SST 305 K on a large domain (6400x400 km)"
end

# ╔═╡ fe212f12-a070-463e-bb85-27c4da4e8665
begin
	fol0Ds = datadir("RCEMIP","small","0D")
	fol1Ds = datadir("RCEMIP","small","1D")
	fol2Ds = datadir("RCEMIP","small","2D")
	md"Folders containing RCEMIP data for SST 305 K on a small domain (96x96 km)"
end

# ╔═╡ 1faec812-15b4-4ce2-8990-ff3643595116
begin
	dsl = NCDataset(joinpath(fol2Dl,"SAM_CRM_RCE_large305_2D_prw.nc"))
	xl  = dsl["x"][1536:2048] / 1000
	yl  = dsl["y"][:] / 1000
	pwl = dsl["prw"][1536:2048,:,end]
	close(dsl)
	
	dsl = NCDataset(joinpath(fol2Dl,"SAM_CRM_RCE_large305_2D_sprw.nc"))
	swv = dsl["sprw"][1536:2048,:,end]
	close(dsl)
	
	csf = pwl ./ swv * 100; csf = csf[:]
	md"Loading 2D humidity data for large-domain RCEMIP simulation ..."
end

# ╔═╡ 4e1e1658-799a-436c-820a-240be56ea3ab
begin
	dss = NCDataset(joinpath(fol2Ds,"SAM_CRM_RCE_small305_2D_prw.nc"))
	xs  = dss["x"][:] / 1000
	ys  = dss["y"][:] / 1000
	pws = dss["prw"][:]
	close(dss)
	md"Loading 2D humidity data for small-domain RCEMIP simulation ..."
end

# ╔═╡ 82a21af2-a3e7-41dc-98fa-da3ff20c5877
begin
	ds3 = NCDataset(joinpath(fol3Dl,"SAM_CRM_RCE_large305_3D_0000540000.nc"))
	l3l = ds3["p"][:]
	rhl = ds3["hur"][1536:2048,:,:]; rhl = reshape(rhl,:,length(l3l))
	close(ds3)
	md"Loading 3D humidity data for large-domain RCEMIP simulation ..."
end

# ╔═╡ 859099a3-e7e9-42aa-8ed7-690c0fdc780a
begin
	d3s = NCDataset(joinpath(fol1Ds,"SAM_CRM_RCE_small305_1D_hur_avg.nc"))
	l3s = d3s["p"][:]
	rhs = d3s["hur_avg"][:,end]
	close(d3s)
	md"Loading 3D humidity data for small-domain RCEMIP simulation ..."
end

# ╔═╡ fd155da0-f454-4c17-94a9-124fcf3a28d0
begin
	arr = [
		[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1],
		[2,2,5,0,0,0,0,0,3,3,3,3,7,0,0,0,0,4,4,6]]
	pplt.close(); f,a = pplt.subplots(arr,aspect=5,axwidth=6,wspace=0,sharey=0)
	
	c = a[1].contourf(xl,yl,pwl[:,:,end]',levels=30:80,extend="both",cmap="Blues")
	a[1].format(
		ylabel="Y / km",ltitle="(a) Large-Domain RCE",
		xlim=(4864,6144),xlocator=(4864:256:6144),
		ylim=(0,256),ylocator=0:128:256
	)
	
	c = a[3].contourf(xs,ys,pws[:,:,end]',levels=30:80,extend="both",cmap="Blues")
	a[3].format(
		xlabel="X / km",ylabel="Y / km",ltitle="(b) Small-Domain RCE",
		xlim=(0,96),xlocator=0:32:96,
		ylim=(0,96),ylocator=0:32:96
	)
	
	a[2].text(0.5,0.5,L"$w_{wtg}$")
	a[2].format(xlocator=[],ylocator=[],ltitle="(c) Dry Regime")
	a[4].text(0.5,0.5,L"$w_{wtg}$")
	a[4].format(xlocator=[],ylocator=[],ltitle="(d) Wet Regime")
	
	ind = sortperm(csf)
	a[5].plot(dropdims(mean(rhl[ind[1:1000],:]*100,dims=1),dims=1),l3l)
	a[5].axis("off")
	a[5].format(ylim=(1010,50),xlim=(0,100),yscale="log")
	
	a[6].plot(dropdims(mean(rhl[ind[(end-1000):end],:]*100,dims=1),dims=1),l3l)
	a[6].axis("off")
	a[6].format(ylim=(1010,50),xlim=(0,100),yscale="log")
	
	a[7].plot(rhs[:],l3s)
	a[7].axis("off")
	a[7].format(ylim=(1010,59),xlim=(0,100),yscale="log")
	
	f.colorbar(c,loc="r",title="Precipitable Water Vapour / mm",locator=30:10:80)
	f.savefig(plotsdir("testRCEMIP.png"),transparent=false,dpi=250)
	load(plotsdir("testRCEMIP.png"))
end

# ╔═╡ Cell order:
# ╟─3f630f34-1d64-11ec-073c-316ff2a18be2
# ╟─5e772289-f828-4876-8456-edc9824716c6
# ╟─5647a8b2-ab77-4cab-b1e0-3d94f0c6d332
# ╟─fe212f12-a070-463e-bb85-27c4da4e8665
# ╟─1faec812-15b4-4ce2-8990-ff3643595116
# ╟─4e1e1658-799a-436c-820a-240be56ea3ab
# ╟─82a21af2-a3e7-41dc-98fa-da3ff20c5877
# ╟─859099a3-e7e9-42aa-8ed7-690c0fdc780a
# ╟─fd155da0-f454-4c17-94a9-124fcf3a28d0
