using DrWatson
@quickactivate "2023GL104350"
using Logging
using Printf

schname = "VDD"
expname = "T1282km300V64"

tprm  = projectdir("exp","tmp.prm")
avec = [1  1   1   1    1   1   1  1 0.5 0]
bvec = [0 0.1 0.2 0.25 0.3 0.4 0.5 1  1  1]
nexp = length(avec)

for iexp in 1 : nexp
    astr  = replace(@sprintf("%04.2f",avec[iexp]),"."=>"d")
    bstr  = replace(@sprintf("%04.2f",bvec[iexp]),"."=>"d")
    conii = "H$(astr)F$(bstr)"
    mkpath(projectdir("exp","prm",schname,expname,conii))
    for imember = 1 : 15
        mstr = @sprintf("%02d",imember)
        oprm = projectdir("run","modifysam","prmtemplates","$schname.prm")
        nprm = projectdir("exp","prm",schname,expname,conii,"member$(mstr).prm")
        open(tprm,"w") do fprm
            open(oprm,"r") do rprm
                s = read(rprm,String)
                s = replace(s,"[expname]" => expname)
                s = replace(s,"[xx]"      => mstr)
                s = replace(s,"[en]"      => "$(imember)")
                s = replace(s,"[aa]"      => "$(avec[iexp])")
                s = replace(s,"[bb]"      => "$(bvec[iexp])")
                write(fprm,s)
            end
        end
        mv(tprm,nprm,force=true)
        @info "Creating new prm file for $schname $expname $conii ensemble member $imember"
    end
end
