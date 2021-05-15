using NCDatasets
using Printf
using Statistics

function createradmean(
    radname::AbstractString,
    raddata::AbstractArray
)

    mkpath(projectdir("exp/rad")); frad = projectdir("exp/rad/$(radname)")
    printrad(frad,raddata)

end

function printrad(frad::AbstractString, raddata::Array{<:Real,2})

    nz = size(raddata,1)

    open(frad,"w") do io
        @printf(io,"p[mb] (dT/dt)rad [K/s]\n")
    end

    open(frad,"a") do io
        for t in [0.0,10000.0]
            @printf(io,"%10.2f, %10d\n",t,nz)
            for iz = 1 : nz
                @printf(
                    io,"%10.4f\t%16.8e\n",
                    raddata[iz,1],raddata[iz,2]
                )
            end
        end
    end

end
