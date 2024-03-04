using DrWatson
@quickactivate "ExploreWTGSpace"
using Logging
using Printf

include(srcdir("sam.jl"))

schname = "DGW"
expname = "P1282km300V64"

if schname == "DGW"
    pwrvec = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
else
    pwrvec = [sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2)]
    pwrvec = vcat(pwrvec/10,1,pwrvec,10,pwrvec*10)
end

tprm  = projectdir("exp","tmp.prm")

for powerii in pwrvec

    pwrname = dampingstrprnt(powerii)
    if schname == "DGW"
          wtgdmp = wtgii; wtgrlx = 1
    else; wtgrlx = wtgii; wtgdmp = 1
    end

    mkpath(projectdir("exp","prm",schname,expname,pwrname))
    for imember = 1 : 15
        mstr = @sprintf("%02d",imember)
        oprm = projectdir("run","modifysam","prmtemplates","$(schname)_$(expname).prm")
        nprm = projectdir("exp","prm",schname,expname,pwrname,"member$(mstr).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"[xx]"  => mstr)
                s = replace(s,"[en]"  => "$(imember)")
                s = replace(s,"[am]"  => @sprintf("%7e",wtgdmp))
                s = replace(s,"[tau]" => @sprintf("%7e",wtgrlx))
                s = replace(s,"e+"    => "e")
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm",schname,expname,pwrname))
        mv(tprm,nprm,force=true)
        @info "Creating new prm file for $schname $expname $pwrname ensemble member $imember"
    end

end
