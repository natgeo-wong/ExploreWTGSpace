using DrWatson
@quickactivate "ExploreWTGSpace"
using Printf

schname = "VDD"
expname = "T1282km300V64"
email   = ""
doBuild = true

tprm  = projectdir("exp","tmp.prm")
avec = [1  1   1   1    1   1   1  1 0.5 0]
bvec = [0 0.1 0.2 0.25 0.3 0.4 0.5 1  1  1]
nexp = length(avec)

mrun = projectdir("run","modifysam","runtemplates","modelrun.sh")
brun = projectdir("run","modifysam","runtemplates","Build.csh")

open(mrun,"r") do frun
    s = read(frun,String)
    for iexp in 1 : nexp
        astr  = replace(@sprintf("%04.2f",avec[iexp]),"."=>"d")
        bstr  = replace(@sprintf("%04.2f",bvec[iexp]),"."=>"d")
        runname = "H$(astr)F$(bstr)"

        for ensembleii in 1 : 15

            mstr = @sprintf("%02d",ensembleii)
            nrun = projectdir("run",schname,expname,runname,"ensemble$(mstr).sh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[email]"       => email)
                sn = replace(sn,"[directory]"   => projectdir())
                sn = replace(sn,"[experiment]"  => expname)
                sn = replace(sn,"[config]"      => runname)
                if ensembleii < 6
                    sn = replace(sn,"[sndname]" => "$(expname)")
                elseif (ensembleii > 5) && (ensembleii < 11)
                    sn = replace(sn,"[sndname]" => "$(expname)_hot")
                else
                    sn = replace(sn,"[sndname]" => "$(expname)_cold")
                end
                sn = replace(sn,"[lsfname]"     => "noforcing")
                sn = replace(sn,"[schname]"     => schname)
                sn = replace(sn,"member[xx]"    => "member$(mstr)")
                write(wrun,sn)
            end

        end

    end
end

if doBuild
    open(brun,"r") do frun
        s = read(frun,String)
        for iexp in 1 : nexp

            astr  = replace(@sprintf("%04.2f",avec[iexp]),"."=>"d")
            bstr  = replace(@sprintf("%04.2f",bvec[iexp]),"."=>"d")
            runname = "H$(astr)F$(bstr)"
            nrun = projectdir("run",schname,expname,runname,"Build.csh")

            open(nrun,"w") do wrun
                sn = replace(s ,"[datadir]" => datadir())
                sn = replace(sn,"[schname]" => schname)
                sn = replace(sn,"[expname]" => expname)
                sn = replace(sn,"[runname]" => runname)
                write(wrun,sn)
            end

        end
    end
end