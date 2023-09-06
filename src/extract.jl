using Dates
using DrWatson
using Logging
using NCDatasets
using Printf
using Statistics

function extractprecip(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    prcp = zeros(nt,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        fnc = datadir(
            "$schname","$expname","$runname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        )
        if isfile(fnc)
            ods = NCDataset(fnc)
            try
                prcp[:,ids] .= ods["PREC"][:]
            catch
                @warn "Unable to extract precipitation data from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    fnc = datadir("precipitation","$(schname)-$(expname)-$(runname).nc")
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"] = nt
    nds.dim["ensemble"] = nmember

    nctime = defVar(nds,"time",Float64,("time",),attrib = Dict(
        "units"     => "days after model-day 0",
        "full_name" => "Day"
    ))

    ncvar = defVar(nds,"precipitation",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "mm day**-1",
        "long_name" => "daily_precipitation_rate",
        "full_name" => "Daily Precipitation Rate"
    ))

    nctime[:] = t
    ncvar[:]  = prcp

    close(nds)

end

function extractprecip(
    schname :: String,
    expname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    prcp = zeros(nt,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        @info "$(Dates.now()) - Opening $(datadir(
            "$schname","$expname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))"
        ods = NCDataset(datadir(
            "$schname","$expname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        ))
        nit = ods.dim["time"]
        prcp[1:nit,ids] .= ods["PREC"][:]
        close(ods)

    end

    fnc = datadir("precipitation","$(schname)-$(expname).nc")
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"] = nt
    nds.dim["ensemble"] = nmember

    nctime = defVar(nds,"time",Float64,("time",),attrib = Dict(
        "units"     => "days after model-day 0",
        "full_name" => "Day"
    ))

    ncvar = defVar(nds,"precipitation",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "mm day**-1",
        "long_name" => "daily_precipitation_rate",
        "full_name" => "Daily Precipitation Rate"
    ))

    nctime[:] = t
    ncvar[:]  = prcp

    close(nds)

end

function extractwwtg(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    wwtg = zeros(64,nt,nmember) * NaN
    z = zeros(64,nmember) * NaN
    p = zeros(64,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    for ids = 1 : nmember

        ensemble = @sprintf("%02d",ids)
        fnc = datadir(
            "$schname","$expname","$runname","OUT_STAT",
            "$(schname)_ExploreWTGSpace-$(expname)-member$(ensemble).nc"
        )
        if isfile(fnc)
            ods = NCDataset(fnc)
            try
                z[:,ids] .= ods["z"][:]
                p[:,ids] .= ods["p"][:]
                wwtg[:,:,ids] .= ods["WWTG"][:]
            catch
                @warn "Unable to extract WWTG data from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    fnc = datadir("wwtg","$(schname)-$(expname)-$(runname).nc")
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"] = nt
    nds.dim["level"] = 64
    nds.dim["ensemble"] = nmember

    nctime = defVar(nds,"time",Float64,("time",),attrib = Dict(
        "units"     => "days after model-day 0",
        "full_name" => "Day"
    ))

    ncz = defVar(nds,"z",Float64,("level",),attrib = Dict(
        "units"     => "m",
        "full_name" => "height"
    ))

    ncp = defVar(nds,"p",Float64,("level","ensemble"),attrib = Dict(
        "units"     => "hPa",
        "full_name" => "pressure_level"
    ))

    ncwwtg = defVar(nds,"wwtg",Float64,("level","time","ensemble"),attrib = Dict(
        "units"     => "m s**-1",
        "long_name" => "weak_temperature_gradient_vertical_velocity",
        "full_name" => "WTG Vertical Velocity"
    ))

    nctime[:] = t
    ncz[:] = dropdims(mean(z,dims=2),dims=2)
    ncp[:] = p
    ncwwtg[:]  = wwtg

    close(nds)

end