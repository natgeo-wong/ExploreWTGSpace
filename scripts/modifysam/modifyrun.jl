using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

include(srcdir("sam.jl"))

schname = "DGW"
expname = "P1282km300V64"
email   = ""
doBuild = true

if checkschemeDGW(schname)
    pwrvec = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
else
    pwrvec = [sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2)]
    pwrvec = vcat(pwrvec/10,1,pwrvec,10,pwrvec*10)
end

mrun = projectdir("run","modifysam","runtemplates","modelrun.sh")
brun = projectdir("run","modifysam","runtemplates","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for pwrii in pwrvec

        pwrname = powername(pwrii,schname)

        for ensembleii in 1 : 15

            mstr = @sprintf("%02d",ensembleii)
            nrun = projectdir("run",schname,expname,pwrname,"ensemble$(mstr).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"   => email)
                sn = replace(sn,"[dirname]" => projectdir())
                sn = replace(sn,"[schname]" => schname)
                sn = replace(sn,"[expname]" => expname)
                sn = replace(sn,"[pwrname]" => pwrname)
                if ensembleii < 6
                    sn = replace(sn,"[sndname]" => "$(expname)")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]" => "$(expname)_warm")
                else
                    sn = replace(sn,"[sndname]" => "$(expname)_cool")
                end
                sn = replace(sn,"[lsfname]" => "noforcing")
                sn = replace(sn,"[memberx]" => "member$(mstr)")
                write(wrun,sn)
            end

        end

    end
end

if doBuild
    open(brun,"r") do frun
        s = read(frun,String)
        for pwrii in pwrvec

            pwrname = powername(pwrii,schname)
            
            nrun = projectdir("run",schname,expname,pwrname,"Build.csh")
            open(nrun,"w") do wrun
                sn = replace(s ,"[datadir]" => datadir())
                sn = replace(sn,"[schname]" => schname)
                sn = replace(sn,"[expname]" => expname)
                sn = replace(sn,"[runname]" => pwrname)
                write(wrun,sn)
            end

        end
    end
end