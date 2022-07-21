using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

expname = "P1282km300V64"
pwrvec1 = [1,2,4,8,16,32,64,128,256,512]
pwrvec2 = [1,2,4,6,8,11,16,23,32,45,64,90,128,256,512]
pwrvec  = pwrvec1

mrun = projectdir("run","modelrun.sh")
brun = projectdir("run","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = "damping$(@sprintf("%04.1f",pwrii))"
        pwrname = replace(pwrname,"."=>"d")

        for ensembleii in 1 : 15

            ensname = "ensemble$(@sprintf("%02d",ensembleii))"
            nrun = projectdir("run",schname,expname,pwrname,"$(ensname).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"=>"")
                sn = replace(s ,"[user]"=>"")
                sn = replace(sn,"[project]"=>"ExploreWTGSpace")
                sn = replace(sn,"[experiment]"=>"$(expname)")
                sn = replace(sn,"[config]"=>"$(pwrname)")
                sn = replace(sn,"[sndname]"=>"$(expname)")
                sn = replace(sn,"[lsfname]"=>"noforcing")
                sn = replace(sn,"[schname]"=>"DGW")
                sn = replace(sn,"member[xx]"=>"member$(@sprintf("%02d",ensembleii))")
                write(wrun,sn)
            end

        end

    end
end

open(brun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = "damping$(@sprintf("%04.1f",pwrii))"
        pwrname = replace(pwrname,"."=>"d")
        nrun = projectdir("run",schname,expname,pwrname,"Build.csh")

        open(nrun,"w") do wrun
            sn = replace(s ,"[user]"=>"")
            sn = replace(sn,"[schname]"=>"DGW")
            sn = replace(sn,"[expname]"=>"$(expname)")
            sn = replace(sn,"[runname]"=>"$(pwrname)")
            write(wrun,sn)
        end

    end
end