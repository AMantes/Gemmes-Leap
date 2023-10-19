library(readxl)
library(xtable)
library(readr)

# Set the working directory with forward slashes
setwd("C:/Users/Achilleas/Documents/GitHub/Gemmes-Leap")
setwd("~/GitHub/Gemmes-Leap/GEMMES")


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



resAlt.d <- as.data.frame(resAlt)

# Create a logical vector to identify annual observations
is_annual <- resAlt.d$time %% 1 == 0  # This checks if time is a whole number (an annual observation)

# Subset the dataframe to keep only the annual observations
annual_resAlt.d <- resAlt.d[is_annual, ]


# Calculate the growth rate of VPIB
#annual_g$g <- c(0, diff(annual_resAlt.d$VPIB) / lag(annual_resAlt.d$VPIB))

# Convert each column to a separate vector
VPIB <- as.vector(annual_resAlt.d$VPIB)
GrowthRateVPIB <- as.vector(annual_resAlt.d$GrowthRateVPIB)

# Create the "t" variable
t <- seq(2007, length.out = length(VPIB))

# Transpose the vectors into a matrix
combined_vector <- rbind(t, VPIB)


#write.csv(resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/resAlt.csv" ,row.names = FALSE)
#write.csv(annual_resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/annual_resAlt.csv" ,row.names = FALSE)
#write.csv(annual_g$g, "C:/Users/Achilleas/Documents/Gemmes/annual_g.csv" ,row.names = FALSE)
#write.csv(transposed_data, "C:/Users/Achilleas/Documents/Gemmes/annual_gt.csv" ,row.names = FALSE)
#write.csv(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = FALSE,col.names = FALSE)
write.table(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = TRUE,col.names = FALSE,dec=".",sep=",")

