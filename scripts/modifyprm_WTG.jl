using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

expname = "S1284km300V64"
tprm  = projectdir("exp","tmp.prm")
plist1 = [
    2,2*sqrt(2.5),2*sqrt(sqrt(2.5))^3,5,5*sqrt(sqrt(2)),
    5*sqrt(2),5*sqrt(sqrt(2))^3,10,10*sqrt(sqrt(2)),10*sqrt(2),
    10*sqrt(sqrt(2))^3,20,20*sqrt(sqrt(2.5)),20*sqrt(2.5),50
]
plist2 = [2,2*sqrt(2.5),5,5*sqrt(2),10,10*sqrt(2),20,20*sqrt(2.5),50]

for powerii in plist1
    conii = "relaxscale$(@sprintf("%04.1f",powerii))"
    conii = replace(conii,"."=>"d")
    mkpath(projectdir("exp","prm","WTG",expname,conii))
    for imember = 1 : 15
        oprm  = projectdir("exp","prm","WTG",expname,"relaxscale02d0","member01.prm")
        nprm  = projectdir("exp","prm","WTG",expname,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,
                    "ttheta_wtg =   2.000000,"=>
                    "ttheta_wtg = $(@sprintf("%10f",powerii)),"
                )
                s = replace(s,"member01"=>"member$(@sprintf("%02d",imember))")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm","WTG",expname,conii))
        mv(tprm,nprm,force=true)
    end
end
