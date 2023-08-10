library(readxl)
library(xtable)

source("SourceCode.R")
source("sourceCodeCalibration.R")
source("utilities.R")

modeleBaseline <- cppMakeMinDist("modelForAchilleas.R", reportVars = 3)

parms<-as.numeric(read.table("parms.csv",sep=";",dec=",")[,1])
names(parms)<-row.names(read.table("parms.csv",sep=";",dec=","))
y0<-as.numeric(read.table("y0.csv",sep=";",dec=",",row.names = 1)[,1])
names(y0)<-row.names(read.table("y0.csv",sep=";",dec=","))

resBaseline <- cppRK4(modeleBaseline,parms=parms,y0=y0)

dataEXOG<-read.csv("defaultEXOG.csv",sep=";",dec=",")

LEAP=read.csv("LEAP_extract.csv",sep=";",dec=",")
dataEXOG$VISOEI_LEAP<-LEAP$VISOEI_LEAP

resAlt <- cppRK4(modeleBaseline,parms=parms,y0=y0,dataExogVar = dataEXOG)

allRes<-list()
allRes[['baseline']]<-resBaseline
allRes[['alternative']]<-resAlt

mymatplotcompare(allRes,"growth(VPIB)","topleft")

# par(mfcol=c(2,4))
# mymatplot(resBaseline,c("VSKNA/VPRODMNA"),"bottomright")
# mymatplot(resBaseline,c("(VYDMNANet-VYDMNAe)/VYDMNAe"),"bottomright")
# mymatplot(resBaseline,c("VSNA/VPIB"),"topright")
# mymatplot(resBaseline,c("PRODMNA/PIB"),"topright")
# mymatplot(resBaseline,c("sigmaMC","sigmaMCC"),"topright")
# mymatplot(resBaseline,c("sigmaMI","sigmaMIC"),"topright")
# mymatplot(resBaseline,c("sigmaMUI","sigmaMUIC"),"topright")
# mymatplot(resBaseline,c("sigmaMX","sigmaMXC"),"topright")
