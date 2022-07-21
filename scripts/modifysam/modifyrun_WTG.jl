using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

schname = "CON"
expname = "P1282km300V64"
tauvec1 = [
    2,2*sqrt(2.5),2*sqrt(sqrt(2.5))^3,5,5*sqrt(sqrt(2)),
    5*sqrt(2),5*sqrt(sqrt(2))^3,10,10*sqrt(sqrt(2)),10*sqrt(2),
    10*sqrt(sqrt(2))^3,20,20*sqrt(sqrt(2.5)),20*sqrt(2.5),50
]
tauvec2 = [2,2*sqrt(2.5),5,5*sqrt(2),10,10*sqrt(2),20,20*sqrt(2.5),50]
tauvec  = tauvec1

mrun = projectdir("run","modelrun.sh")
brun = projectdir("run","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for tauii in tauvec

        tauname = "relaxscale$(@sprintf("%04.1f",tauii))"
        tauname = replace(tauname,"."=>"d")

        for ensembleii in 1 : 15

            ensname = "ensemble$(@sprintf("%02d",ensembleii))"
            nrun = projectdir("run",schname,expname,tauname,"$(ensname).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"=>"")
                sn = replace(sn,"[user]"=>"")
                sn = replace(sn,"[project]"=>"ExploreWTGSpace")
                sn = replace(sn,"[experiment]"=>"$(expname)")
                sn = replace(sn,"[config]"=>"$(tauname)")
                if ensembleii < 6
                    sn = replace(sn,"[sndname]"=>"$(expname)")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]"=>"$(expname)_hot")
                else
                    sn = replace(sn,"[sndname]"=>"$(expname)_cld")
                end
                sn = replace(sn,"[lsfname]"=>"noforcing")
                sn = replace(sn,"[schname]"=>"$(schname)")
                sn = replace(sn,"member[xx]"=>"member$(@sprintf("%02d",ensembleii))")
                write(wrun,sn)
            end

        end

    end
end

open(brun,"r") do frun
    s = read(frun,String)
    for tauii in tauvec

        tauname = "relaxscale$(@sprintf("%04.1f",tauii))"
        tauname = replace(tauname,"."=>"d")
        nrun = projectdir("run",schname,expname,tauname,"Build.csh")

        open(nrun,"w") do wrun
            sn = replace(s ,"[user]"=>"")
            sn = replace(sn,"[schname]"=>"$(schname)")
            sn = replace(sn,"[expname]"=>"$(expname)")
            sn = replace(sn,"[runname]"=>"$(tauname)")
            write(wrun,sn)
        end

    end
end