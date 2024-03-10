using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

expname = "P1282km300V64"
tprm  = projectdir("exp","tmp.prm")
mkpath(projectdir("exp","prm","RCE",expname))

for imember = 1 : 10
    mstr = @sprintf("%02d",imember)
    oprm  = projectdir("run","modifysam","prmtemplates","RCE_$(expname).prm")
    nprm  = projectdir("exp","prm","RCE",expname,"member$(mstr).prm")
    open(tprm,"w") do fprm
        open(oprm,"r") do rprm
            s = read(rprm,String)
            s = replace(s,"[xx]" => mstr)
            s = replace(s,"[en]" => "$(imember)")
            write(fprm,s)
        end
    end
    mv(tprm,nprm,force=true)
end