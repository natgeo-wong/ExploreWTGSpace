using DrWatson
@quickactivate "ExploreWTGSpace"

include(srcdir("sam.jl"))
include(srcdir("extract.jl"))

schname = "DGW"
expname = "P1282km300V64"

if schname == "DGW"
    pwrvec = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
else
    pwrvec = [sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2)]
    pwrvec = vcat(TGR/10,1,TGR,10,TGR*10)
end

for pwr in pwrvec
    pwrname = powername(pwr,schname)
    extractprecip(schname,expname,pwrname,nt=6000,tperday=24)
    extractwwtg(schname,expname,pwrname,nt=6000,tperday=24)
    extractgms(schname,expname,pwrname,nt=6000,tperday=24)
end