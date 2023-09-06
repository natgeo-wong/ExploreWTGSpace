using DrWatson
@quickactivate "2023GL104350"
using Printf

include(srcdir("samsnd.jl"))
include(srcdir("sam.jl"))

expii = "P1282km300V64"

z,p,_,q,t,_ = readsnd(projectdir("exp","snd",expii))
nz = length(z)
snddata = zeros(nz,6)
snddata[:,1] .= z
snddata[:,2] .= p
snddata[:,4] .= q

t_h = t .+ 0.05; snddata[:,3] .= t_h .* (1000 ./p).^(287/1004)
createsndmean("$(expii)_hot",snddata;psfc=1009.32)

t_h = t .- 0.05; snddata[:,3] .= t_h .* (1000 ./p).^(287/1004)
createsndmean("$(expii)_cld",snddata;psfc=1009.32)