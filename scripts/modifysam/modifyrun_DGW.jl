using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

include(srcdir("sam.jl"))

expname = "P1282km300V64"
pwrvec1 = [0,0.2,0.5,1,2,5,10,20,50,100,200,500]
pwrvec2 = [0,0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
pwrvec  = pwrvec1

mrun = projectdir("run","modifysam","runtemplates","modelrun.sh")
brun = projectdir("run","modifysam","runtemplates","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = dampingstrprnt(pwrii)

        for ensembleii in 1 : 15

            mstr = @sprintf("%02d",ensembleii)
            nrun = projectdir("run","DGW",expname,pwrname,"ensemble$(mstr).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"       => )
                sn = replace(sn,"[directory]"   => projectdir())
                sn = replace(sn,"[experiment]"  => expname)
                sn = replace(sn,"[config]"      => pwrname)
                if ensembleii < 6
                    sn = replace(sn,"[sndname]" => "$(expname)")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]" => "$(expname)_hot")
                else
                    sn = replace(sn,"[sndname]" => "$(expname)_cld")
                end
                sn = replace(sn,"[lsfname]"     => "noforcing")
                sn = replace(sn,"[schname]"     => "DGW")
                sn = replace(sn,"member[xx]"    => "member$(mstr)")
                write(wrun,sn)
            end

        end

    end
end

open(brun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = dampingstrprnt(pwrii)
        nrun = projectdir("run","DGW",expname,pwrname,"Build.csh")

        open(nrun,"w") do wrun
            sn = replace(s ,"[datadir]" => datadir())
            sn = replace(sn,"[schname]" => "DGW")
            sn = replace(sn,"[expname]" => expname)
            sn = replace(sn,"[runname]" => pwrname)
            write(wrun,sn)
        end

    end
end