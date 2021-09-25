using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

exp = "D0641km300"
tprm = projectdir("exp","tmp.prm")

for powerii = [1,2,4,8,16,32,64,128,256,512]
    conii = "damping$(@sprintf("%03d",powerii))"
    mkpath(projectdir("exp","prm",exp,conii))
    for imember = 1 : 10
        prm  = projectdir("exp","prm",exp,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s = read(oprm,String)
                s = replace(s,"tabs_s     = 301.70,"=>"tabs_s     = 300.,")
                write(fprm,s)
            end
        end
        mv(tprm,prm,force=true)
    end
end
