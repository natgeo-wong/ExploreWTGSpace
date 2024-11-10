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

    fol = datadir("precipitation","$(schname)","$(expname)")
    fnc = joinpath(fol,"$(runname).nc")
    if !isdir(fol); mkpath(fol) end
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
    ncvar[:,:]  = prcp

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

    fol = datadir("precipitation","$(schname)")
    fnc = joinpath(fol,"$(expname).nc")
    if !isdir(fol); mkpath(fol) end
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
    ncvar[:,:]  = prcp

    close(nds)

end

function extractwwtg(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
    nmode   :: Int = 10
)

    wwtg = zeros(64,nt,nmember) * NaN
    mode = zeros(nt,nmode,nmember)
    ztrop = zeros(nt,nmember)
    ptrop = zeros(nt,nmember)
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
                wwtg[:,:,ids] .= ods["WWTG"][:,:]
            catch
                @warn "Unable to extract WWTG data from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    for ids = 1 : nmember, it = 1 : nt

        iwwtg = @views wwtg[:,it,ids]
        if !isnan(sum(iwwtg))
            itrop = findlast(iwwtg!=0) + 1
            ztrop[it,ids] = z[itrop,ids]
            ptrop[it,ids] = p[itrop,ids]
            iz = @views z[1:itrop,ids]
            for imode = 1 : nmode
                iiwwtg = iwwtg .* sin.(iz./ztrop[it,ids]*pi*imode)
                mode[it,imode,ids] = 2/ztrop[it,ids]*trapz(vcat(0,iz),vcat(0,iiwwtg))
            end
        else
            mode[it,:,ids] .= NaN
        end

    end

    fol = datadir("wwtg","$(schname)","$(expname)")
    fnc = joinpath(fol,"$(runname).nc")
    if !isdir(fol); mkpath(fol) end
    if isfile(fnc); rm(fnc,force=true) end

    nds = NCDataset(fnc,"c",attrib = Dict(
        "Conventions" => "CF-1.6",
        "history"     => "Created on $(Dates.now()) using NCDatasets.jl",
        "comments"    => "Creating NetCDF files in the same format that data is saved on the Climate Data Store"
    ))

    nds.dim["time"]  = nt
    nds.dim["level"] = 64
    nds.dim["mode"]  = nmode
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

    ncmode = defVar(nds,"wₙ",Float64,("time","mode","ensemble"),attrib = Dict(
        "units"     => "m s**-1",
        "long_name" => "weak_temperature_gradient_vertical_mode_coefficient",
        "full_name" => "WTG Vertical Mode Coefficient"
    ))

    ncztrop = defVar(nds,"ztrop",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "m",
        "long_name" => "tropopause_height",
        "full_name" => "Tropopause Height"
    ))

    ncptrop = defVar(nds,"ptrop",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "hPa",
        "long_name" => "tropopause_pressure",
        "full_name" => "Tropopause Pressure"
    ))

    nctime[:] = t
    ncz[:] = dropdims(mean(z,dims=2),dims=2)
    ncp[:] = p
    ncwwtg[:,:,:] = wwtg
    ncmode[:,:,:] = mode
    ncztrop[:] = ztrop
    ncptrop[:] = ptrop

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
                w[:,:,ids] .= ods["WWTG"][:,:]
                dse[:,:,ids]  .= ods["DSE"][:,:]
                mse[:,:,ids]  .= ods["MSE"][:,:]
            catch
                @warn "Unable to extract WWTG data from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    fol = datadir("gms","$(schname)","$(expname)")
    fnc = joinpath(fol,"$(runname).nc")
    if !isdir(fol); mkpath(fol) end
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
    ncdse[:,:,:] = dse
    ncmse[:,:,:] = mse

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

    ncwdse[:,:] = wdse
    ncwmse[:,:] = wmse
    ncgms[:,:]  = wmse ./ wdse

    close(nds)

end

function extractmoisturemode(
    schname :: String,
    expname :: String,
    runname :: String;
    nt      :: Int = 2000,
    tperday :: Int = 1,
    nmember :: Int = 15,
)
    
    z = zeros(64,nmember) * NaN
    t = 0 : nt
    t = collect(t[1:(end-1)] .+ t[2:end]) / 2
    t = t / tperday

    pre = zeros(64,nt,nmember) * NaN
    tbs = zeros(64,nt,nmember) * NaN
    qbs = zeros(64,nt,nmember) * NaN

    int_tbs = zeros(nt,nmember) * NaN
    int_qbs = zeros(nt,nmember) * NaN

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
                pre[:,:,ids] .= ods["PRES"][:,:] * 100
                tbs[:,:,ids] .= ods["TBIAS"][:,:]
                qbs[:,:,ids] .= ods["QBIAS"][:,:]
            catch
                @warn "Unable to extract statistical data for moisture mode calculations from $(fnc)"
            end
            close(ods)
        else
            @warn "No file exists at $(fnc), please run this configuration again ..."
        end

    end

    ipp = zeros(66); ipp[1] = 100923; iipp  = @view ipp[2:(end-1)]
    itmp = zeros(66); iitmp = @view itmp[2:(end-1)]

    for ids = 1 : nmember, it = 1 : nt

        iipp  .= pre[:,it,ids]
        iitmp .= tbs[:,it,ids]; int_tbs[it,ids] = -trapz(ipp,itmp) * 1.0035
        iitmp .= qbs[:,it,ids]; int_qbs[it,ids] = -trapz(ipp,itmp) * 2.5009

    end

    fol = datadir("moisturemode","$(schname)","$(expname)")
    fnc = joinpath(fol,"$(runname).nc")
    if !isdir(fol); mkpath(fol) end
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

    ncmse = defVar(nds,"hp",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "J kg**-1",
        "long_name" => "anomalous_column_integrated_moist_static_energy",
        "full_name" => "Anomalous Column Integrated Moist Static Energy"
    ))

    nctbs = defVar(nds,"Tp",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "J kg**-1",
        "long_name" => "anomalous_column_integrated_temperature",
        "full_name" => "Anomalous Column Integrated Temperature"
    ))

    ncqbs = defVar(nds,"qp",Float64,("time","ensemble"),attrib = Dict(
        "units"     => "J kg**-1",
        "long_name" => "anomalous_column_integrated_water_vapour",
        "full_name" => "Anomalous Column Integrated Water Vapour"
    ))

    nctime[:] = t
    ncz[:] = dropdims(mean(z,dims=2),dims=2)
    nctbs[:,:] = int_tbs
    ncqbs[:,:] = int_qbs
    ncmse[:,:] = int_tbs .+ int_qbs

    close(nds)

end