library(readxl)
library(xtable)
library(readr)
library(openxlsx)
library(reshape2)



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
LEAP=read.csv("inv_costs_clean.csv",sep=",",dec=".") #cleaned_foo


LEAP$diff <- LEAP$netzero-LEAP$reference
LEAP_inv <- LEAP$diff[26: length(LEAP$diff)]
plot(LEAP$diff,type="l")

#LEAP = t(LEAP)
#LEAP = as.data.frame(LEAP)
#LEAP=as.numeric(LEAP)
#LEAP=as.data.frame(LEAP)
#LEAP = LEAP[-c(1),] # non numeric
#LEAP = as.numeric(LEAP)
#LEAP = as.data.frame(LEAP)
#LEAP = LEAP[-c(1:25),] #
#LEAP = as.numeric(LEAP)
#LEAP = as.data.frame(LEAP)

dataEXOG$VISOEI_LEAP<-LEAP_inv

####################
#LEAP=read.csv("LEAP_extract.csv",sep=",",dec=".") #cleaned_foo
#LEAP=as.data.frame(t(LEAP))
#dataEXOG$VISOEI_LEAP<-LEAP$Total
#dataEXOG$VISOEI_LEAP<-LEAP$VISOEI_LEAP

resAlt <- cppRK4(modeleBaseline,parms=parms,y0=y0,dataExogVar = dataEXOG)


allRes<-list()
allRes[['baseline']]<-resBaseline
allRes[['alternative']]<-resAlt

#mymatplotcompare(allRes,"growth(VPIB)","topleft")


resAlt.d <- as.data.frame(resAlt)

# Create a logical vector to identify annual observations
is_annual <- resAlt.d$time %% 1 == 0  # This checks if time is a whole number (an annual observation)

# Subset the dataframe to keep only the annual observations
annual_resAlt.d <- resAlt.d[is_annual, ]


# Calculate the growth rate of VPIB
#annual_g$g <- c(0, diff(annual_resAlt.d$VPIB) / lag(annual_resAlt.d$VPIB))

# Convert each column to a separate vector
VPIB <- as.vector(annual_resAlt.d$VPIB)
VPIB<-c(1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1044.96,1050.41,1103.54,1137.37,1.015131*1.2386*VPIB/1000,2789.402604,2878.663487,2970.780719,3065.845702,3163.952764,3265.199253,3369.685629,3477.515569,3588.796067,3703.637541)*1000

GrowthRateVPIB <- as.vector(annual_resAlt.d$GrowthRateVPIB)

# Create the "t" variable
t <- seq(2007, length.out = length(VPIB))

# Transpose the vectors into a matrix
combined_vector <- rbind(t, VPIB)


write.csv(resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/resAlt.csv" ,row.names = FALSE)
#write.csv(annual_resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/annual_resAlt.csv" ,row.names = FALSE)
#write.csv(annual_g$g, "C:/Users/Achilleas/Documents/Gemmes/annual_g.csv" ,row.names = FALSE)
#write.csv(transposed_data, "C:/Users/Achilleas/Documents/Gemmes/annual_gt.csv" ,row.names = FALSE)
#write.csv(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = FALSE,col.names = FALSE)
#write.table(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = TRUE,col.names = FALSE,dec=".",sep=",")
combined_vector_df <- as.data.frame(combined_vector)

xlsx_file_path <- "C:/Users/Achilleas/Documents/Gemmes/combined_vector.xlsx"
reshaped_df <- melt(combined_vector, id.vars = NULL)
combined_vector_df1 <- data.frame(t = combined_vector[1, ], VPIB = combined_vector[2, ])
#write.xlsx(combined_vector_df1, xlsx_file_path, rowNames = TRUE,sheet = "Sheet")
# Specify the path where you want to save the XLSX file
xlsx_file_path <- "C:/Users/Achilleas/Documents/Gemmes/combined_vector.xlsx"
# Create a workbook and add a data frame to it with the sheet name "Sheet"
wb <- createWorkbook()
addWorksheet(wb, "Sheet")
writeData(wb, sheet = 1, x = combined_vector_df1)
# Save the workbook to the specified path
saveWorkbook(wb, xlsx_file_path,overwrite = TRUE)

# par(mfcol=c(2,4))
# mymatplot(resBaseline,c("VSKNA/VPRODMNA"),"bottomright")
# mymatplot(resBaseline,c("(VYDMNANet-VYDMNAe)/VYDMNAe"),"bottomright")
# mymatplot(resBaseline,c("VSNA/VPIB"),"topright")
# mymatplot(resBaseline,c("PRODMNA/PIB"),"topright")
# mymatplot(resBaseline,c("sigmaMC","sigmaMCC"),"topright")
# mymatplot(resBaseline,c("sigmaMI","sigmaMIC"),"topright")
# mymatplot(resBaseline,c("sigmaMUI","sigmaMUIC"),"topright")
# mymatplot(resBaseline,c("sigmaMX","sigmaMXC"),"topright")
