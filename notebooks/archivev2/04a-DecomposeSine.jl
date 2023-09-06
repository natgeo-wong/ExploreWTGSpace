### A Pluto.jl notebook ###
# v0.19.16

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
	using PlutoUI
	using Printf
	using StatsBase
	
	using ImageShow, PNGFiles
	using PyCall, LaTeXStrings
	pplt = pyimport("proplot")
	
	include(srcdir("sam.jl"))
	
md"Loading modules for the ExploreWTGSpace project..."
end

# ╔═╡ e78a75c2-590f-11eb-1144-9127b0309135
md"
# 4a. Fourier Decomposition of the WTG

In this notebook, we explore our attempts to merge the WTG temperature relaxation schemes with those of the Damped Gravity Wave.
"

# ╔═╡ b6892634-9199-11eb-38d5-8dda8da16ed7
md"
### A. Comparing the DGW and WTG Schemes

We recall the momentum damping equation (assuming $a_m$ is constant with height):

$$\frac{\partial^2\omega'}{\partial p^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'$$

Let us rewrite this equation into vertical coordinates:

$$\frac{\partial^2\omega'}{\partial z^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v'
\cdot \frac{\partial^2 p}{\partial z^2} \approx \frac{k^2}{a_m} \frac{R_d}{p}T_v' (\rho g)^2 \approx \frac{k^2}{a_m} \frac{\rho^2g^2}{\rho T} T_v' = \frac{k^2}{a_m} \frac{pg^2}{RT^2} T_v'$$

From the fact that $\omega = -w\rho g$, we obtain

$$\therefore (\rho w')_{zz} \approx -\frac{k^2}{a_m} \frac{\rho g}{T} T_v'$$

We assume that the half-sine and full-sine modes dominate both $\rho w'$, $\theta'$ and $T_v'$. Furthermore, we assume that $\theta'$ and $T_v'$ are roughly equal in magnitude, and that hydrostatic balance is obeyed, such that:

$$\rho w' = \rho_0 e^{-z/H} (w_1 \sin(\pi z / z_t) + w_2 \sin(2\pi z / z_t))$$
$$\theta' \approx T_v' = T_1 \sin(\pi z / z_t) + T_2 \sin(2\pi z / z_t)$$

The second derivative of $\rho w'$ by height is given by

$$(\rho w')_{zz}(1) \approx w_1\rho_0 ((1/H^2-\pi^2/z_t^2)e^{-z/H}\sin(\pi z/z_t) - 2\pi/Hz_t e^{-z/H}\cos(\pi z/z_t))$$
$$(\rho w')_{zz}(2) \approx w_2\rho_0 ((1/H^2-4\pi^2/z_t^2)e^{-z/H}\sin(2\pi z/z_t) - 4\pi/Hz_t e^{-z/H}\cos(2\pi z/z_t))$$

We assume that the $\cos$ components do not factor into our analysis, which leaves us with the following

$$(\rho w')_{zz} = \rho (w_1(1/H^2-\pi^2/z_t^2)\sin(\pi z/z_t) + w_2(1/H^2-4\pi^2/z_t^2)\sin(2\pi z/z_t))$$

So we have that

$$\rho w_1(1/H^2-\pi^2/z_t^2)\sin(\pi z/z_t) = -\frac{k^2}{a_m} \frac{\rho g}{T} T_1 \sin(\pi z / z_t)$$
$$\rho w_2(1/H^2-4\pi^2/z_t^2)\sin(2\pi z/z_t) = -\frac{k^2}{a_m} \frac{\rho g}{T} T_2 \sin(2\pi z / z_t)$$

Such that

$$w_1 = \frac{k^2}{a_m} \frac{g}{T(\pi^2/z_t^2-1/H^2)} T_1$$
$$w_2 = \frac{k^2}{a_m} \frac{g}{T(4\pi^2/z_t^2-1/H^2)} T_2$$

Comparing this to the Weak Temperature Gradient scheme, given that we assume that $\partial\theta/\partial z = 0.001$ K m$^{-1}$, gives us

$$w_1 = \frac{k^2}{a_m} \frac{g}{T(\pi^2/z_t^2-1/H^2)} T_1 = \frac{1000T_1}{\tau_1}$$

Thus, we see that

$$\tau_1 = \frac{1000T(\pi^2/z_t^2-1/H^2)}{gk^2}a_m$$

$$\tau_2 = \frac{1000T(4\pi^2/z_t^2-1/H^2)}{gk^2}a_m$$

So we see that in order to approximate the Damped Gravity Wave impacts, we must adjust such that $\tau_2 \approx 4\tau_1$. And if we want to add additional wave components, we will add an $n^2$ to each relevant relaxation factor.
"

# ╔═╡ 67da2e88-2139-4fb9-b79f-3e0e2f91f5df
md"
### B. Varying $\tau_h$ and $\tau_f$ in the WTG Schemes

We now run a series of simulations that vary either $\tau_f$ or $\tau_h$, while keeping the other fixed at a given value of $\tau$, which in this case is 15 hr.  We vary the strength of $\tau_f$ and $\tau_h$ such that the ratio of $\tau/\tau_f$ and $\tau/\tau_h$ are varied from 0 to 1. When $\tau/\tau_f = 0$ this means that the full-sine component of $\theta'$ has zero impact on the final $w'$ (i.e. the resultant $w'$ is fully dominated by the half-sine component of $\theta'$), and vice-versa.
"

# ╔═╡ 55230f4a-7661-11eb-1c37-8b022b95e08e
# ╠═╡ show_logs = false
begin
	pplt.close()
	fts,ats = pplt.subplots(
		[1,1,1,2,3,3,3,4],aspect=1.5,axwidth=2,wspace=[0,0,0,3,0,0,0],sharex=0
	)

	for imem = 1 : 10
		fnc = outstatname("RCE","S1284km300V64","",false,true,imem)
		if isfile(fnc)
			pr = retrievevar_fnc("PREC",fnc)./24
			ats[1].plot([0,1],[1,1]*mean(pr[(end-499):end]),c="grey",lw=1)
			ats[2].plot([0,1],[1,1]*mean(pr[(end-499):end]),c="grey",lw=1)
			ats[3].plot([0,1],[1,1]*mean(pr[(end-499):end]),c="grey",lw=1)
			ats[4].plot([0,1],[1,1]*mean(pr[(end-499):end]),c="grey",lw=1)
		end
	end

	for ic in 0 : 0.1 : 0.5
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")

		imem = 0
		while imem < 15; imem += 1
			fnc = outstatname("VDD",expname,"damping010",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				pr = mean(retrievevar_fnc("PREC",fnc))./24

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[1].plot(ic,pr,marker=".",c="k",ms=4)
				else
					ats[1].scatter(ic,pr,c=clr,alpha=0.1,s=50)
				end
			end
		end
		
	end

	for ic in 1
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")

		imem = 0
		while imem < 15; imem += 1
			fnc = outstatname("VDD",expname,"damping010",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				pr = mean(retrievevar_fnc("PREC",fnc))./24

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[1].plot(0.6,pr,marker=".",c="k",ms=4)
				else
					ats[1].scatter(0.6,pr,c=clr,alpha=0.1,s=50)
				end
			end
		end
		
	end

	for ic in 0 : 0.1 : 1
		expname = "H$(@sprintf("%4.2f",ic))F1d00"
        expname = replace(expname,"."=>"d")

		imem = 0
		while imem < 15; imem += 1
			fnc = outstatname("VDD",expname,"damping010",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				pr = mean(retrievevar_fnc("PREC",fnc))./24

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[2].plot(ic,pr,marker=".",c="k",ms=4)
				else
					ats[2].scatter(ic,pr,c=clr,alpha=0.1,s=50)
				end
			end
		end
		
	end

	

	for ic in 0 : 0.1 : 0.5
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")

		imem = 0
		while imem < 15; imem += 1
			fnc = outstatname("VDT",expname,"relaxscale10d0",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				pr = mean(retrievevar_fnc("PREC",fnc))./24

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[3].plot(ic,pr,marker=".",c="k",ms=4)
				else
					ats[3].scatter(ic,pr,c=clr,alpha=0.1,s=50)
				end
			end
		end
		
	end

	for ic in 1
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")

		imem = 0
		while imem < 15; imem += 1
			fnc = outstatname("VDT",expname,"relaxscale10d0",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				pr = mean(retrievevar_fnc("PREC",fnc))./24

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[3].plot(0.6,pr,marker=".",c="k",ms=4)
				else
					ats[3].scatter(0.6,pr,c=clr,alpha=0.1,s=50)
				end
			end
		end
		
	end

	for ic in 0 : 0.1 : 1
		expname = "H$(@sprintf("%4.2f",ic))F1d00"
        expname = replace(expname,"."=>"d")

		imem = 0
		while imem < 15; imem += 1
			fnc = outstatname("VDT",expname,"relaxscale10d0",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				pr = mean(retrievevar_fnc("PREC",fnc))./24

				if (imem >= 6) && (imem <= 10)
					clr = "yellow7"
				else
					clr = "blue5"
				end

				if imem <= 5
					ats[4].plot(ic,pr,marker=".",c="k",ms=4)
				else
					ats[4].scatter(ic,pr,c=clr,alpha=0.1,s=50)
				end
			end
		end

	end

	for iats = 1 : 2 : 4
		ats[iats].format(
			xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
		)
		ats[iats+1].format(
			yloc="right",xlim=(1,0),
			xlabel=L"\tau/\tau_h",xlocator=0:0.5:1,
		)
	end

	ats[1].format(
		xlabel=L"c_2",title=L"$a_m=10$ day$^{-1}$",
		ultitle=L"(a) $c_1 = 1$",
		xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
	)
	ats[2].format(
		xlim=(1,0),urtitle=L"(b) $c_2 = 1$",
		xlabel=L"c_1",xlocator=0:0.5:1,
	)

	ats[3].format(
		xlabel=L"c_2",title=L"$\tau=10$ hr",
		ultitle=L"(a) $c_1 = 1$",
		xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
	)
	ats[4].format(
		xlim=(1,0),urtitle=L"(b) $c_2 = 1$",
		xlabel=L"c_1",xlocator=0:0.5:1,
	)
	ats[2].format(ylim=(0,1),ylabel=L"Rain Rate / mm hr$^{-1}$")
	
	fts.savefig(plotsdir("04a-decomposesine.png"),transparent=false,dpi=400)
	load(plotsdir("04a-decomposesine.png"))
end

# ╔═╡ 9cf4fa56-91a8-11eb-2710-955eefd10142
md"
We see the following:
* As $\tau/\tau_f$ decreases, the dry and wet regime-states visible in the DGW simulations start to become more prominent.
* As $\tau/\tau_h$ decreases, the final solution becomes more similar to RCE

As shown by our calculations above, $\tau/\tau_f \approx 1/4$ in the DGW simulations, and there we see that the separation between the dry and wet regime-states are noticeable.
"

# ╔═╡ 9cfdf73b-0f59-4b7b-9a30-cba33d9c343e
begin
	pplt.close()
	fpr,apr = pplt.subplots(
		[[1,1,1,2,3,3,3,4,5,5,5,6],[7,7,7,8,9,9,9,10,11,11,11,12]],
		aspect=1,axwidth=1.2,wspace=[0,0,0,2,0,0,0,2,0,0,0],sharex=0
	)

	wwtg = zeros(64,7,2); wwtg_f = zeros(64,3,2)
	hsin = zeros(64,7,2); hsin_f = zeros(64,3,2)
	fsin = zeros(64,7,2); fsin_f = zeros(64,3,2)
	ppre = zeros(64,7,2); ppre_f = zeros(64,3,2)

	for ic in 0 : 0.1 : 0.5
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")
		icon = Int(ic*10+1)

		nmem = 0
		for imem = 6 : 10
			fnc = outstatname("VDD",expname,"damping050",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				nmem += 1
				_,p,_ = retrievedims_fnc(fnc)
				ppre[:,icon,1] += p
				wwtg[:,icon,1] += dropdims(mean(retrievevar_fnc("WWTG",fnc),dims=2),dims=2)
				hsin[:,icon,1] += dropdims(mean(retrievevar_fnc("WWTGHSIN",fnc),dims=2),dims=2)
				fsin[:,icon,1] += dropdims(mean(retrievevar_fnc("WWTGFSIN",fnc),dims=2),dims=2)
			end
		end
		ppre[:,icon,1] = ppre[:,icon,1] / nmem
		wwtg[:,icon,1] = wwtg[:,icon,1] / nmem
		hsin[:,icon,1] = hsin[:,icon,1] / nmem
		fsin[:,icon,1] = fsin[:,icon,1] / nmem

		nmem = 0
		for imem = 11 : 15
			fnc = outstatname("VDD",expname,"damping050",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				nmem += 1
				_,p,_ = retrievedims_fnc(fnc)
				ppre[:,icon,2] += p
				wwtg[:,icon,2] += dropdims(mean(retrievevar_fnc("WWTG",fnc),dims=2),dims=2)
				hsin[:,icon,2] += dropdims(mean(retrievevar_fnc("WWTGHSIN",fnc),dims=2),dims=2)
				fsin[:,icon,2] += dropdims(mean(retrievevar_fnc("WWTGFSIN",fnc),dims=2),dims=2)
			end
		end
		ppre[:,icon,2] = ppre[:,icon,2] / nmem
		wwtg[:,icon,2] = wwtg[:,icon,2] / nmem
		hsin[:,icon,2] = hsin[:,icon,2] / nmem
		fsin[:,icon,2] = fsin[:,icon,2] / nmem
		
	end

	for ic in 0 : 0.5 : 1
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")
		icon = Int(ic*2+1)

		nmem = 0
		for imem = 6 : 10
			fnc = outstatname("VDD",expname,"damping050",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				nmem += 1
				_,p,_ = retrievedims_fnc(fnc)
				ppre_f[:,icon,1] += p
				wwtg_f[:,icon,1] += dropdims(mean(retrievevar_fnc("WWTG",fnc),dims=2),dims=2)
				hsin_f[:,icon,1] += dropdims(mean(retrievevar_fnc("WWTGHSIN",fnc),dims=2),dims=2)
				fsin_f[:,icon,1] += dropdims(mean(retrievevar_fnc("WWTGFSIN",fnc),dims=2),dims=2)
			end
		end
		ppre_f[:,icon,1] = ppre_f[:,icon,1] / nmem
		wwtg_f[:,icon,1] = wwtg_f[:,icon,1] / nmem
		hsin_f[:,icon,1] = hsin_f[:,icon,1] / nmem
		fsin_f[:,icon,1] = fsin_f[:,icon,1] / nmem

		nmem = 0
		for imem = 11 : 15
			fnc = outstatname("VDD",expname,"damping050",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				nmem += 1
				_,p,_ = retrievedims_fnc(fnc)
				ppre_f[:,icon,2] += p
				wwtg_f[:,icon,2] += dropdims(mean(retrievevar_fnc("WWTG",fnc),dims=2),dims=2)
				hsin_f[:,icon,2] += dropdims(mean(retrievevar_fnc("WWTGHSIN",fnc),dims=2),dims=2)
				fsin_f[:,icon,2] += dropdims(mean(retrievevar_fnc("WWTGFSIN",fnc),dims=2),dims=2)
			end
		end
		ppre_f[:,icon,2] = ppre_f[:,icon,2] / nmem
		wwtg_f[:,icon,2] = wwtg_f[:,icon,2] / nmem
		hsin_f[:,icon,2] = hsin_f[:,icon,2] / nmem
		fsin_f[:,icon,2] = fsin_f[:,icon,2] / nmem
		
	end

	for ic in 1
		expname = "H1d00F$(@sprintf("%4.2f",ic))"
        expname = replace(expname,"."=>"d")

		nmem = 0
		for imem = 6 : 10
			fnc = outstatname("VDD",expname,"damping050",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				nmem += 1
				_,p,_ = retrievedims_fnc(fnc)
				ppre[:,7,1] += p
				wwtg[:,7,1] += dropdims(mean(retrievevar_fnc("WWTG",fnc),dims=2),dims=2)
				hsin[:,7,1] += dropdims(mean(retrievevar_fnc("WWTGHSIN",fnc),dims=2),dims=2)
				fsin[:,7,1] += dropdims(mean(retrievevar_fnc("WWTGFSIN",fnc),dims=2),dims=2)
			end
		end
		ppre[:,7,1] = ppre[:,7,1] / nmem
		wwtg[:,7,1] = wwtg[:,7,1] / nmem
		hsin[:,7,1] = hsin[:,7,1] / nmem
		fsin[:,7,1] = fsin[:,7,1] / nmem

		nmem = 0
		for imem = 11 : 15
			fnc = outstatname("VDD",expname,"damping050",false,true,imem)
			fnc = replace(fnc,"-$expname-"=>"-S1284km300V64-")
			if isfile(fnc)
				nmem += 1
				_,p,_ = retrievedims_fnc(fnc)
				ppre[:,7,2] += p
				wwtg[:,7,2] += dropdims(mean(retrievevar_fnc("WWTG",fnc),dims=2),dims=2)
				hsin[:,7,2] += dropdims(mean(retrievevar_fnc("WWTGHSIN",fnc),dims=2),dims=2)
				fsin[:,7,2] += dropdims(mean(retrievevar_fnc("WWTGFSIN",fnc),dims=2),dims=2)
			end
		end
		ppre[:,7,2] = ppre[:,7,2] / nmem
		wwtg[:,7,2] = wwtg[:,7,2] / nmem
		hsin[:,7,2] = hsin[:,7,2] / nmem
		fsin[:,7,2] = fsin[:,7,2] / nmem
		
	end

	wwtg = wwtg * 1000; wwtg_f = wwtg_f * 1000
	hsin = hsin * 1000; hsin_f = hsin_f * 1000
	fsin = fsin * 1000; fsin_f = fsin_f * 1000

	wwtg_f[:,3,:] = wwtg[:,7,:]
	hsin_f[:,3,:] = hsin[:,7,:]
	fsin_f[:,3,:] = fsin[:,7,:]

	pdry = dropdims(mean(ppre[:,:,1],dims=2),dims=2)
	pwet = dropdims(mean(ppre[:,:,2],dims=2),dims=2)

	lvls = vcat(-50,-20,-10,-5,-2,-1,-0.5,0.5,1,2,5,10,20,50)/10
	
	apr[1].contourf(0:0.1:0.6,pdry,wwtg[:,:,1],cmap="drywet",levels=lvls,extend="both")
	apr[2].contourf(0:0.5:1,pdry,wwtg_f[:,:,1],cmap="drywet",levels=lvls,extend="both")
	apr[3].contourf(0:0.1:0.6,pdry,hsin[:,:,1],cmap="drywet",levels=lvls,extend="both")
	apr[4].contourf(0:0.5:1,pdry,hsin_f[:,:,1],cmap="drywet",levels=lvls,extend="both")
	apr[5].contourf(0:0.1:0.6,pdry,fsin[:,:,1],cmap="drywet",levels=lvls,extend="both")
	apr[6].contourf(0:0.5:1,pdry,fsin_f[:,:,1],cmap="drywet",levels=lvls,extend="both")

	c = apr[7].contourf( 0:0.1:0.6,pwet,wwtg[:,:,2],cmap="drywet",levels=lvls,extend="both")
	apr[8].contourf( 0:0.5:1,pwet,wwtg_f[:,:,2],cmap="drywet",levels=lvls,extend="both")
	apr[9].contourf( 0:0.1:0.6,pwet,hsin[:,:,2],cmap="drywet",levels=lvls,extend="both")
	apr[10].contourf( 0:0.5:1,pwet,hsin_f[:,:,2],cmap="drywet",levels=lvls,extend="both")
	apr[11].contourf(0:0.1:0.6,pwet,fsin[:,:,2],cmap="drywet",levels=lvls,extend="both")
	apr[12].contourf(0:0.5:1,pwet,fsin_f[:,:,2],cmap="drywet",levels=lvls,extend="both")

	for iats = 1 : 2 : 12
		apr[iats].format(
			suptitle=L"$\tau=15$ hr",
			xlocator=0:0.1:0.5,xminorticks=0:0.02:0.5,xlim=(0,0.6)
		)
		apr[iats+1].format(yloc="right",xlim=(1,0),xlocator=0:0.5:1,)
	end

	apr[1].format(ylim=(1000,25),yscale="log",ylabel="Pressure / hPa")
	apr[7].format(ylim=(1000,25),yscale="log",ylabel="Pressure / hPa",xlabel=L"\tau/\tau_f")
	apr[8].format(xlabel=L"\tau/\tau_h")
	apr[9].format(xlabel=L"\tau/\tau_f")
	apr[10].format(xlabel=L"\tau/\tau_h")
	apr[11].format(xlabel=L"\tau/\tau_f")
	apr[12].format(xlabel=L"\tau/\tau_h")

	apr[1].format(ltitle=L"(a) Total $w_{WTG}$",leftlabels=["Dry Regime","Wet Regime"])
	apr[3].format(ltitle=L"(b) Half-Sine $w_{WTG}$")
	apr[5].format(ltitle=L"(c) Full-Sine $w_{WTG}$")

	fpr.colorbar(c,length=0.75,label=L"$10^{-3}$ m s$^{-1}$")
	fpr.savefig(plotsdir("04a-decomposesine-vertprofiles.png"),transparent=false,dpi=400)
	load(plotsdir("04a-decomposesine-vertprofiles.png"))
end

# ╔═╡ Cell order:
# ╟─e78a75c2-590f-11eb-1144-9127b0309135
# ╟─681658b0-5914-11eb-0d65-bbace277d145
# ╟─6dce35fc-5914-11eb-0ce2-0d4e164e1898
# ╟─b6892634-9199-11eb-38d5-8dda8da16ed7
# ╟─67da2e88-2139-4fb9-b79f-3e0e2f91f5df
# ╠═55230f4a-7661-11eb-1c37-8b022b95e08e
# ╟─9cf4fa56-91a8-11eb-2710-955eefd10142
# ╟─9cfdf73b-0f59-4b7b-9a30-cba33d9c343e
