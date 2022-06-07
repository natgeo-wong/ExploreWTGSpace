using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

expname = "S1284km300V64"
tprm  = projectdir("exp","tmp.prm")
plist1 = [1,2,4,8,16,32,64,128,256,512]
plist2 = [1,2,4,6,8,10,12,16,20,24,28,32,40,48,56,64,80,96,112,128,160,192,224,256,512]
plist3 = [1,2,4,6,8,11,16,23,32,45,64,90,128,256,512]

for powerii in plist1
    conii = "damping$(@sprintf("%03d",powerii))"
    mkpath(projectdir("exp","prm","DGW",expname,conii))
    for imember = 1 : 15
        oprm  = projectdir("exp","prm","DGW",expname,"damping001","member01.prm")
        nprm  = projectdir("exp","prm","DGW",expname,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"am_wtg  = 1., "=>"am_wtg  = $(powerii)., ")
                s = replace(s,"member01"=>"member$(@sprintf("%02d",imember))")
                s = replace(s,"nensemble     = 1,"=>"nensemble     = $(imember),")
                write(fprm,s)
            end
        end
        mkpath(projectdir("exp","prm","DGW",expname,conii))
        mv(tprm,nprm,force=true)
    end
end
