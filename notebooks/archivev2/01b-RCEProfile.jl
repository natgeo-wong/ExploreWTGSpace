### A Pluto.jl notebook ###
# v0.16.4

using Markdown
using InteractiveUtils

# ╔═╡ df810659-6173-483f-912c-523b022d641e
begin
	using Pkg; Pkg.activate()
	using DrWatson
	md"Using DrWatson to ensure reproducibility between different machines ..."
end

# ╔═╡ 46faa412-5ade-11eb-3c37-23a7e59037a0
begin
	@quickactivate "ExploreWTGSpace"
	using DelimitedFiles
	using Statistics
	using Printf
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
	md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ 9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
md"
# 1b. Plotting the RCE Profiles

This notebook extracts the sounding profile data for different RCE experiments so that we can compare and contrast the RCE profiles for different configurations.
"

# ╔═╡ b2670c08-81e5-11eb-324e-2b923b289a04
md"
### A. Loading the Vertical Profile Data ...
"

# ╔═╡ bcc23ba4-c40e-4966-94d9-326feddf0b12
begin
	expvec = [
		"P0641km295V64","P0641km300V64","P0641km305V64",
		"P1281km300V64",
		"P1282km295V64","P1282km300V64","P1282km305V64",
		"P1284km300V64",
		# "P1284km300","S1284km300",
		# "S0641km300",
		# "T0641km300",
	]
	colvec = [
		"b","g","r",
		"green8",
		"blue6","green6","red6",
		"green3",
		"green1","blue3","blue2","blue1"]
	nexp = length(expvec)
end

# ╔═╡ c13062a4-373b-489d-a619-e7b33f059ba7
begin
	pplt.close()
	f3D,a3D = pplt.subplots(ncols=6,aspect=0.4,axwidth=1.2,sharex=0)
	lgd = Dict("ncol"=>1,"frame"=>false)
	
	for ie in 1 : nexp
		
		iexp = expvec[ie]
		imem = 0
		
		while imem < 100; imem += 1
			
			fnc = outstatname("Control",expvec[ie],false,true,imem)
	
			if isfile(fnc)
				
				_,p,t = retrievedims(fnc); t = t .- 80
				clci = mean(retrievevar("CLD",fnc)[:,(end-99):end],dims=2)*100
				tabi = mean(retrievevar("TABS",fnc)[:,(end-99):end],dims=2)
				qvi  = mean(retrievevar("QV",fnc)[:,(end-99):end],dims=2) / 10
				rhi  = calcrh(qvi,tabi,p)
				
				nets = mean(retrievevar("RADQRSW",fnc)[:,(end-99):end],dims=2)
				netl = mean(retrievevar("RADQRLW",fnc)[:,(end-99):end],dims=2)
				neth = mean(retrievevar("RADQR",fnc)[:,(end-99):end],dims=2)
				
				a3D[1].plot(dropdims(clci,dims=2),p,lw=1,c=colvec[ie])
				a3D[2].plot(dropdims(tabi,dims=2),p,lw=1,c=colvec[ie])
				a3D[3].plot(dropdims(rhi,dims=2),p,lw=1,c=colvec[ie])
				a3D[4].plot(dropdims(nets,dims=2),p,lw=1,c=colvec[ie])
				a3D[5].plot(dropdims(netl,dims=2),p,lw=1,c=colvec[ie])
				
				if imem == 1
				a3D[6].plot(dropdims(neth,dims=2),p,lw=1,c=colvec[ie],label=expvec[ie],legend="r",legend_kw=lgd,)
				else
					a3D[6].plot(dropdims(neth,dims=2),p,lw=1,c=colvec[ie])
				end
				
			end
			
		end
	
	end
	
	a3D[1].format(xlim=(0,25),xlabel="Cloud Cover / %")
	a3D[2].format(xlim=(180,320),xlabel="Temperature / K")
	a3D[3].format(xlim=(0,110),xlabel="Relative Humidity / %")
	a3D[4].format(xlim=(-2.5,2.5),xlabel=L"SW Heating / K day$^{-1}$")
	a3D[5].format(xlim=(-2.5,2.5),xlabel=L"LW Heating / K day$^{-1}$")
	a3D[6].format(xlim=(-2.5,2.5),xlabel=L"Net Heating / K day$^{-1}$")
	
	for ax in a3D
		ax.format(ylim=(1000,20),yscale="log")
	end
	
	f3D.savefig(plotsdir(
		"01b-rceprofile-3D.png"),
		transparent=false,dpi=200
	)
	load(plotsdir("01b-rceprofile-3D.png"))
end

# ╔═╡ 967cf880-bf6c-4f3b-80ef-c477bc4cdbaf


# ╔═╡ 0a102885-12c6-433e-9999-db457517bd0b


# ╔═╡ Cell order:
# ╟─9dd4cd7e-5adb-11eb-2735-a7a4a2bb23b1
# ╟─df810659-6173-483f-912c-523b022d641e
# ╟─46faa412-5ade-11eb-3c37-23a7e59037a0
# ╟─b2670c08-81e5-11eb-324e-2b923b289a04
# ╠═bcc23ba4-c40e-4966-94d9-326feddf0b12
# ╠═c13062a4-373b-489d-a619-e7b33f059ba7
# ╠═967cf880-bf6c-4f3b-80ef-c477bc4cdbaf
# ╠═0a102885-12c6-433e-9999-db457517bd0b
