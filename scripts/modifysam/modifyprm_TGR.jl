using DrWatson
@quickactivate "2023GL104350"
using Logging
using Printf

include(srcdir("sam.jl"))

expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")
pvec  = [sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2)]
pvec1 = vcat(pwrvec,10,pwrvec*10)
pvec2 = vcat(pwrvec/10,1,pwrvec,10,pwrvec*10)
pvec  = pvec1

for powerii in pvec
    conii = relaxscalestrprnt(powerii)
    mkpath(projectdir("exp","prm","TGR",expname,conii))
    for imember = 1 : 15
        mstr = @sprintf("%02d",imember)
        oprm = projectdir("run","modifysam","prmtemplates","TGR_$(expname).prm")
        nprm = projectdir("exp","prm","TGR",expname,conii,"member$(mstr).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"[xx]" => mstr)
                s = replace(s,"[en]" => "$(imember)")
                s = replace(s,"[am]" => @sprintf("%7e",1))
                if !iszero(powerii)
                    s = replace(s,"[bool]" => "true")
                    s = replace(s,"[tau]" => @sprintf("%7e",powerii))
                else
                    s = replace(s,"[bool]" => "false")
                    s = replace(s,"[tau]" => @sprintf("%7e",1))
                end
                s = replace(s,"e+" => "e")
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm","TGR",expname,conii))
        mv(tprm,nprm,force=true)
        @info "Creating new prm file for TGR $expname $conii ensemble member $imember"
    end
end
