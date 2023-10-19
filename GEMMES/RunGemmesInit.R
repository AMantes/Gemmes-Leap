library(readxl)
library(xtable)


# Set the working directory with forward slashes
setwd("C:/Users/Achilleas/Documents/GitHub/Gemmes-Leap")
setwd("~/GitHub/Gemmes-Leap/GEMMES")


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
#####################
LEAP=read.csv("cleaned_foo2.csv",sep=",",dec=".") #cleaned_foo
LEAP = t(LEAP)
#LEAP = as.data.frame(LEAP)
#LEAP=as.numeric(LEAP)
#LEAP=as.data.frame(LEAP)
LEAP = LEAP[-c(1),] # non numeric
#LEAP = as.numeric(LEAP)
LEAP = as.data.frame(LEAP)
LEAP = LEAP[-c(1:25),] #
LEAP = as.numeric(LEAP)
LEAP = as.data.frame(LEAP)

dataEXOG$VISOEI_LEAP<-LEAP$LEAP

####################
#LEAP=read.csv("LEAP_extract.csv",sep=",",dec=".") #cleaned_foo
#LEAP=as.data.frame(t(LEAP))
#dataEXOG$VISOEI_LEAP<-LEAP$Total
#dataEXOG$VISOEI_LEAP<-LEAP$VISOEI_LEAP

resAlt <- cppRK4(modeleBaseline,parms=parms,y0=y0,dataExogVar = dataEXOG)


allRes<-list()
allRes[['baseline']]<-resBaseline
allRes[['alternative']]<-resAlt

mymatplotcompare(allRes,"growth(VPIB)","topleft")


resAlt.d <- as.data.frame(resAlt)

# Create a logical vector to identify annual observations
is_annual <- resAlt.d$time %% 1 == 0  # This checks if time is a whole number (an annual observation)

# Subset the dataframe to keep only the annual observations
annual_resAlt.d <- resAlt.d[is_annual, ]


# Calculate the growth rate of VPIB
annual_g$g <- c(0, diff(annual_resAlt.d$VPIB) / lag(annual_resAlt.d$VPIB))

# Convert each column to a separate vector
VPIB <- as.vector(annual_resAlt.d$VPIB)
GrowthRateVPIB <- as.vector(annual_resAlt.d$GrowthRateVPIB)

# Create the "t" variable
t <- seq(2007, length.out = length(VPIB))

# Transpose the vectors into a matrix
combined_vector <- rbind(t, VPIB)


write.csv(resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/resAlt.csv" ,row.names = FALSE)
write.csv(annual_resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/annual_resAlt.csv" ,row.names = FALSE)
write.csv(annual_g$g, "C:/Users/Achilleas/Documents/Gemmes/annual_g.csv" ,row.names = FALSE)
write.csv(transposed_data, "C:/Users/Achilleas/Documents/Gemmes/annual_gt.csv" ,row.names = FALSE)
write.csv(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = FALSE,col.names = FALSE)
write.table(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = TRUE,col.names = FALSE,dec=".",sep=",")

# par(mfcol=c(2,4))
# mymatplot(resBaseline,c("VSKNA/VPRODMNA"),"bottomright")
# mymatplot(resBaseline,c("(VYDMNANet-VYDMNAe)/VYDMNAe"),"bottomright")
# mymatplot(resBaseline,c("VSNA/VPIB"),"topright")
# mymatplot(resBaseline,c("PRODMNA/PIB"),"topright")
# mymatplot(resBaseline,c("sigmaMC","sigmaMCC"),"topright")
# mymatplot(resBaseline,c("sigmaMI","sigmaMIC"),"topright")
# mymatplot(resBaseline,c("sigmaMUI","sigmaMUIC"),"topright")
# mymatplot(resBaseline,c("sigmaMX","sigmaMXC"),"topright")
