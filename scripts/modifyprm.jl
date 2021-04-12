using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

exp = "P064km295d0"
tprm = projectdir("exp","tmp.prm")

for powerii = 0 : 9
    conii = "damping$(@sprintf("%03d",2^powerii))"
    for imember = 1 : 10
        prm = projectdir("exp","prm",exp,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s = read(oprm,String)
                s = replace(s,"tabs_s     = 305.00,"=>"tabs_s     = 295.00,")
                write(fprm,s)
            end
        end
        mv(tprm,prm,force=true)
    end
end
