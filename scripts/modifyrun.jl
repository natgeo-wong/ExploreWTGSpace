using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

exp = "P1282km300V64"
plist1 = [1,2,4,8,16,32,64,128,256,512]

for powerii in plist1
    conii = "damping$(@sprintf("%03d",powerii))"
    trun  = projectdir("run",exp,conii,"ensemblexx.sh")
    open(trun,"r") do frun
        s = read(frun,String)
        for imember = 1 : 10
            nrun = projectdir("run",exp,conii,"ensemble$(@sprintf("%02d",imember)).sh")
            open(nrun,"w") do wrun
                sn = replace(s,"=member00"=>"=member$(@sprintf("%02d",imember))")
                write(wrun,sn)
            end
        end
    end
end
