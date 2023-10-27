library(readxl)
library(xtable)
library(readr)
library(openxlsx)


# Set the working directory with forward slashes
setwd("C:/Users/Achilleas/Documents/GitHub/Gemmes-Leap")
setwd("~/GitHub/Gemmes-Leap/GEMMES")


dataEXOG<-read.csv("defaultEXOG.csv",sep=";",dec=",")
#####################
#LEAP=read.csv("cleaned_foo2.csv",sep=",",dec=".") #cleaned_foo
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
LEAP=read.csv("inv_costs_clean.csv",sep=",",dec=".") #cleaned_foo
LEAP$diff <- LEAP$netzero-LEAP$reference
LEAP_inv <- LEAP$diff[26: length(LEAP$diff)]
plot(LEAP$diff,type="l")

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



resAlt.d <- as.data.frame(resAlt)

# Get the current time in "HHMM" format
current_time <- format(Sys.time(), "%H%M")

# Create a filename with the current time and a .csv extension
filename <- paste("resAlt.d", current_time, ".csv", sep = "")

# Export the dataframe as a CSV
write.csv(resAlt.d, file = filename, row.names = FALSE)




# Create a logical vector to identify annual observations
is_annual <- resAlt.d$time %% 1 == 0  # This checks if time is a whole number (an annual observation)

# Subset the dataframe to keep only the annual observations
annual_resAlt.d <- resAlt.d[is_annual, ]


# Calculate the growth rate of VPIB
#annual_g$g <- c(0, diff(annual_resAlt.d$VPIB) / lag(annual_resAlt.d$VPIB))

# Convert each column to a separate vector
VPIB <- as.vector(annual_resAlt.d$VPIB)
VPIB<-c(1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1001.45,1044.96,1050.41,1103.54,1137.37,1.015131*1.2386*VPIB/1000,2789.402604,2878.663487,2970.780719,3065.845702,3163.952764,3265.199253,3369.685629,3477.515569,3588.796067,3703.637541)*1000
#GrowthRateVPIB <- as.vector(annual_resAlt.d$GrowthRateVPIB)

# Create the "t" variable
t <- seq(2004, length.out = length(VPIB))

# Transpose the vectors into a matrix
combined_vector <- rbind(t, VPIB)


#write.csv(resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/resAlt.csv" ,row.names = FALSE)
#write.csv(annual_resAlt.d, "C:/Users/Achilleas/Documents/Gemmes/annual_resAlt.csv" ,row.names = FALSE)
#write.csv(annual_g$g, "C:/Users/Achilleas/Documents/Gemmes/annual_g.csv" ,row.names = FALSE)
#write.csv(transposed_data, "C:/Users/Achilleas/Documents/Gemmes/annual_gt.csv" ,row.names = FALSE)
#write.csv(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = FALSE,col.names = FALSE)
#write.table(combined_vector, "C:/Users/Achilleas/Documents/Gemmes/combined_vector.csv" ,row.names = TRUE,col.names = FALSE,dec=".",sep=",")


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

#plot(t,VPIB,type = "l")


