using Glob
using NCDatasets
using NumericalIntegration
using Printf
using Statistics

function outstatname(
    expname :: AbstractString,
    config  :: AbstractString,
    iscontrol  :: Bool = true,
    isDGW      :: Bool = false,
    isWTG      :: Bool = false,
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer = 0
)

    if isensemble
    	  expmem = "$(expname)-member$(@sprintf("%02d",member))"
    else; expmem = expname
    end

    prefix = "RCE"
    if isDGW
        prefix = "DGW"
    elseif isWTG
        prefix = "WTG"
    end

    if iscontrol
        fnc = datadir(joinpath(
    		prefix,expname,"OUT_STAT",
    		"$(prefix)_ExploreWTGSpace-$(expmem).nc"
    	))
    elseif ismpi
    	fnc = datadir(joinpath(
    		prefix,"Ensemble",expname,config,"OUT_STAT",
    		"$(prefix)_ExploreWTGSpace-$(expmem).nc"
    	))
    else
    	fnc = datadir(joinpath(
    		prefix,experiment,config,"OUT_STAT",
    		"$(prefix)_ExploreWTGSpace-$(expmem).nc"
    	))
    end

    return fnc

end

function retrievedims(
    expname :: AbstractString,
    config  :: AbstractString = "";
    iscontrol  :: Bool = true,
    isDGW      :: Bool = false,
    isWTG      :: Bool = false,
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds = NCDataset(outstatname(
        expname,config,
        iscontrol,isDGW,isWTG,
        ismpi,isensemble,member
    ))
    z = ds["z"][:]
    p = ds["p"][:]
    t = ds["time"][:]
    close(ds)

    return z,p,t

end

function retrievedims(fnc::AbstractString)

    rce = NCDataset(fnc)
    z = rce["z"][:]
    p = rce["p"][:]
    t = rce["time"][:]
    close(rce)

    return z,p,t

end

function retrievevar(
    varname :: AbstractString,
    expname :: AbstractString,
    config  :: AbstractString;
    iscontrol  :: Bool = true,
    isDGW      :: Bool = false,
    isWTG      :: Bool = false,
    ismpi      :: Bool = false,
    isensemble :: Bool = false,
    member     :: Integer=0
)

    ds = NCDataset(outstatname(
        expname,config,
        iscontrol,isDGW,isWTG,
        ismpi,isensemble,member
    ))
    var = ds[varname][:]
    close(ds)

    return var

end

function retrievevar(variable::AbstractString, fnc::AbstractString)

    ds = NCDataset(fnc)
    var = ds[variable][:]
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
		swp[it] = integrate(pvec,QVsat) / 9.81 / 1000
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
