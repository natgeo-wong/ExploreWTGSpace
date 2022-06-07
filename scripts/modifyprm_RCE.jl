using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")

for imember = 1 : 10
    oprm  = projectdir("exp","prm","RCE",expname,"member01.prm")
    nprm  = projectdir("exp","prm","RCE",expname,"member$(@sprintf("%02d",imember)).prm")
    open(tprm,"w") do fprm
        open(oprm,"r") do rprm
            s = read(rprm,String)
            s = replace(s,"member01"=>"member$(@sprintf("%02d",imember))")
            s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
            write(fprm,s)
        end
    end
    mv(tprm,nprm,force=true)
end
