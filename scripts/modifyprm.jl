using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

exp = "S1282km300V64"
tprm = projectdir("exp","tmp.prm")
plist1 = [1,2,4,8,16,32,64,128,256,512]
plist2 = [1,2,4,6,8,10,12,16,20,24,28,32,40,48,56,64,80,96,112,128,160,192,224,256,512]
plist3 = [1,2,4,6,8,11,16,23,32,45,64,90,128,256,512]

for powerii in plist1
    conii = "damping$(@sprintf("%03d",powerii))"
    mkpath(projectdir("exp","prm",exp,conii))
    for imember = 1 : 10
        prm  = projectdir("exp","prm",exp,conii,"member$(@sprintf("%02d",imember)).prm")
        open(tprm,"w") do fprm
            open(prm,"r") do oprm
                s = read(oprm,String)
                s = replace(s,"dx = 1000.,"=>"dx = 2000.,")
                s = replace(s,"dy = 1000.,"=>"dy = 2000.,")
                write(fprm,s)
            end
        end
        mv(tprm,prm,force=true)
    end
end
