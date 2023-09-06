using DrWatson
@quickactivate "2023GL104350"

include(srcdir("sam.jl"))
include(srcdir("extract.jl"))

schname = "DGW"
expname = "P1282km300V64"

DGW1 = [0.2,0.5,1,2,5,10,20,50,100,200,500]
DGW2 = [0.02,0.05,0.1,0.2,0.5,1,2,5,10,20,50,100,200,500]
TGR  = [sqrt(2),2,2*sqrt(2.5),5,5*sqrt(2)]
TGR1 = vcat(TGR,10,TGR*10)
TGR2 = vcat(TGR/10,1,TGR,10,TGR*10)

pwrvec  = DGW1

for pwr in pwrvec
    extractprecip(schname,expname,dampingstrprnt(pwr),nt=6000,tperday=24)
end