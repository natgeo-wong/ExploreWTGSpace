using NCDatasets
using Printf
using Statistics
using Trapz

dampingstrprnt(am::Real)     = replace(   "damping$(@sprintf("%06.2f",am ))","."=>"d")
relaxscalestrprnt(tau::Real) = replace("relaxscale$(@sprintf("%06.3f",tau))","."=>"d")

function powername(
    power  :: Real,
    scheme :: AbstractString
)

    if (scheme == "DGW") || (scheme == "LGW")
        return dampingstrprnt(power)
    else
        return relaxscalestrprnt(power)
    end

end

function outstatname(
    scheme  :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString,
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer = 0
)

    if isensemble
    	  expmem = "$(expname)-member$(@sprintf("%02d",member))"
    else; expmem = expname
    end

    if scheme == "RCE"
          iscontrol = true
    else; iscontrol = false
    end

    if iscontrol
        fnc = datadir(joinpath(
    		scheme,expname,"OUT_STAT",
    		"$(scheme)_ExploreWTGSpace-$(expmem).nc"
    	))
    elseif ismpi
    	fnc = datadir(joinpath(
    		scheme,"Ensemble",expname,config,"OUT_STAT",
    		"$(scheme)_ExploreWTGSpace-$(expmem).nc"
    	))
    else
    	fnc = datadir(joinpath(
    		scheme,expname,config,"OUT_STAT",
    		"$(scheme)_ExploreWTGSpace-$(expmem).nc"
    	))
    end

    return fnc

end

function retrievedims(
    scheme  :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString = "";
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds = NCDataset(outstatname(
        scheme,expname,config,
        ismpi,isensemble,member
    ))
    z = ds["z"][:]
    p = ds["p"][:]
    t = ds["time"][:]
    close(ds)

    return z,p,t

end

function retrievedims_fnc(fnc::AbstractString)

    rce = NCDataset(fnc)
    z = rce["z"][:]
    p = rce["p"][:]
    t = rce["time"][:]
    close(rce)

    return z,p,t

end

function retrieve2Dvar(
    varname :: AbstractString,
    scheme  :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString = "";
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds = NCDataset(outstatname(
        scheme,expname,config,
        ismpi,isensemble,member
    ))
    var = ds[varname][:]
    close(ds)

    return var

end

function retrieve2Dvar_fnc(varname::AbstractString, fnc::AbstractString)

    ds = NCDataset(fnc)
    var = ds[varname][:]
    close(ds)

    return var

end

function retrieve3Dvar(
    varname :: AbstractString,
    scheme  :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString = "";
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds = NCDataset(outstatname(
        scheme,expname,config,
        ismpi,isensemble,member
    ))
    var = ds[varname][:,:]
    close(ds)

    return var

end

function retrieve3Dvar_fnc(varname::AbstractString, fnc::AbstractString)

    ds = NCDataset(fnc)
    var = ds[varname][:,:]
    close(ds)

    return var

end

function calcrh(QV,TAIR,P)

    RH = zeros(size(QV)); np = size(RH,1); nt = size(RH,2)

    for it = 1 : nt, ip = 1 : np
    	RH[ip,it] = QV[ip,it] / tair2qsat(TAIR[ip,it],P[ip]*100)
    end

    return RH

end

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

function calcswp(RH,QV,P)

	pvec = vcat(0,reverse(P)) * 100; nt = size(RH,2)
	QVsat = zeros(length(pvec))
	swp = zeros(nt)

    for it = 1 : nt
		QVsat[2:end] .= reverse(QV[:,it]) ./ reverse(RH[:,it])
		swp[it] = trapz(pvec,QVsat) / 9.81 / 1000
    end

	return swp

end

function t2d(t::Vector{<:Real}, days::Integer)

    tstep = round(Integer,(length(t)-1)/(t[end]-t[1]))
    t = mod.(t[(end-tstep+1):end],1); tmin = argmin(t)
    tshift = tstep-tmin+1; t = circshift(t,tshift)
    t = vcat(t[end]-1,t,t[1]+1)
    beg = days*tstep - 1

    return t*tstep,tstep,tshift,beg

end
