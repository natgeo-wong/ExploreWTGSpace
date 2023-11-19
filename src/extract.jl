using Dates
using DrWatson
using Logging
using NCDatasets
using Printf
using Statistics
using Trapz

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

function extractgms(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)

    dse = zeros(64,nt,nmember) * NaN
    mse = zeros(64,nt,nmember) * NaN
    wdse = zeros(nt,nmember) * NaN
    wmse = zeros(nt,nmember) * NaN
    w = zeros(64,nt,nmember) * NaN
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
                w[:,:,ids] .= ods["WWTG"][:]
                dse[:,:,ids]  .= ods["DSE"][:]
                mse[:,:,ids]  .= ods["MSE"][:]
            catch
                @warn "Unable to extract WWTG data from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    mkpath(datadir("gms"))
    fnc = datadir("gms","$(schname)-$(expname)-$(runname).nc")
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

    ncdse = defVar(nds,"dse",Float64,("level","time","ensemble"),attrib = Dict(
        "units"     => "J kg**-1",
        "long_name" => "dry_static_energy",
        "full_name" => "Dry Static Energy"
    ))

    ncmse = defVar(nds,"mse",Float64,("level","time","ensemble"),attrib = Dict(
        "units"     => "J kg**-1",
        "long_name" => "moist_static_energy",
        "full_name" => "Moist Static Energy"
    ))

    ncgms = defVar(nds,"gms",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "",
        "long_name" => "gross_moist_stability",
        "full_name" => "Gross Moist Stability"
    ))

    ncwdse = defVar(nds,"wdse",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "J m**2 kg**-1 s**-1",
        "long_name" => "column_integrated_dry_static_energy_export",
        "full_name" => "Column Integrated Dry Static Energy_export"
    ))

    ncwmse = defVar(nds,"wmse",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "J m**2 kg**-1 s**-1",
        "long_name" => "column_integrated_moist_static_energy_export",
        "full_name" => "Column Integrated Moist Static Energy_export"
    ))

    nctime[:] = t
    ncz[:] = dropdims(mean(z,dims=2),dims=2)
    ncp[:] = p
    ncdse[:] = dse
    ncmse[:] = mse

    for imem = 1 : 15

        iz = @view z[:,imem]
        iw = @view w[:,:,imem]
        imse = @view mse[:,:,imem]
        idse = @view dse[:,:,imem]

        dz = (iz[2:end] .- iz[1:(end-1)])
        iz = (iz[2:end] .+ iz[1:(end-1)]) / 2
        iw = (iw[2:end,:] .+ iw[1:(end-1),:]) / 2
        wdsdz = (idse[2:end,:] .- idse[1:(end-1),:]) ./ dz .* iw
        wdhdz = (imse[2:end,:] .- imse[1:(end-1),:]) ./ dz .* iw

        for it = 1 : nt
            wmse[it,imem] = trapz(iz,@views wdhdz[:,it])
            wdse[it,imem] = trapz(iz,@views wdsdz[:,it])
        end

    end

    ncwdse[:] = wdse
    ncwmse[:] = wmse
    ncgms[:]  = wmse ./ wdse

    close(nds)

end