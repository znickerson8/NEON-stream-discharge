######################################################################################################################## 
#' @title Stage-Discharge Rating Curve Controls Script - D02 - LEWI

#' @author Bobby Hensley \email{hensley@battelleecology.org} \cr 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr
#' Nick Harrison \email{nharrison@battelleecology.org} \cr

#' @description This script generates the controls, uncertainties, and priors associated with the creation of a stage-
#' discharge rating curve for Lewis Run for water years 2012-2019.

#' @return This script produces three .csv files:
#' 'geo_controlInfo_in' contains information on the number of controls and their activations
#' 'geo_controlType_in' Defines the control type and reports parameters and uncertainties for each control
#' 'geo_priorParameters_in' reports the priors calculated in this script

#' @references 
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

# changelog and author contributions / copyrights
#   Kaelin Cawley and Nick Harrison (2019)
#     Generic script created
#   Bobby Hensley (2/3/2020)
#     Modified for LEWI 2019 survey
######################################################################################################################## 

#This reads in data using the API and pulls zip files from the ECS buckets
#load packages
library(neonUtilities)
library(plotly)

siteID <- "LEWI"
domainID <- "D02"
streamMorphoDPID <- "DP4.00131.001"
filepath <- "N:/Science/AQU/Controls/D02_LEWI_20190205_new"
URIpath <- paste(filepath,"filesToStack00131","stackedFiles",sep = "/")

# #Download data from API and store somewhere
dataFromAPI <- neonUtilities::zipsByProduct(streamMorphoDPID,siteID,package="expanded",check.size=FALSE,savepath = filepath)
neonUtilities::stackByTable(filepath=paste(filepath,"filesToStack00131",sep = "/"), folder = TRUE)
neonUtilities::zipsByURI(filepath=URIpath, savepath = URIpath, pick.files=FALSE, unzip = TRUE, check.size = FALSE)

#Read in downloaded data
surveyPtsDF <- read.table(paste0(URIpath,"/NEON_D02_LEWI_GEOMORPH_20190205_L0_VE/LEWI_surveyPts_20190205.CSV"),
                          sep = ",",
                          header = TRUE,
                          stringsAsFactors = FALSE)

#The end date of the geomorphology survey (YYYYMMDD)
surveyDate<-'20190205'

#The date when this survey applies to the gauging record
surveyActiveDate <- "2012-01-01" #1/1/2012 is used for the first survey for a site out of convenience

#Stipulate 4-digit site code, underscore, and survey year (ex: HOPB_2017)
surveyID <- "LEWI_2019"

#Creates dataframe of all points associated with transect DSC1
names(surveyPtsDF) <- c("name","latitude","Longitude","northing","easting","elevation","mapCode","E","N","H")
dischargePointsXS1<-subset(surveyPtsDF,mapCode=="Transect_DSC")
dischargePointsXS1<-dischargePointsXS1[order(dischargePointsXS1$N),]
rownames(dischargePointsXS1)<-seq(length=nrow(dischargePointsXS1))

#Sets plot1 settings.  
xAxisTitle1<-list(title="Easting (m)",zeroline=FALSE)
yAxisTitle1<-list(title="Northing (m)",zeroline=FALSE)
font<-list(size=12,color='black')

#Plot the cross section by easting and northing data for a sanity check
plot_ly(data=dischargePointsXS1,x=~E, y=~N, name='Easting vs Northing', type='scatter', mode='markers', text=~name)%>%
  layout(title = siteID, xaxis=xAxisTitle1, yaxis=yAxisTitle1)

#Manually select NorthStart and EastStart coordinates
dischargeXS1NorthStart<-dischargePointsXS1$N[dischargePointsXS1$name=="DSC_LB_PIN"]
dischargeXS1EastStart<-dischargePointsXS1$E[dischargePointsXS1$name=="DSC_LB_PIN"]

#Assigns a raw Distance value to each point relative to the NorthStart and EastStart coordinates.
for(i in 1:(length(dischargePointsXS1$name))){
  dischargeXS1PointN<-dischargePointsXS1$N[i]
  dischargeXS1PointE<-dischargePointsXS1$E[i]
  dischargePointsXS1$DistanceRaw[i]<-sqrt(((dischargeXS1PointN-dischargeXS1NorthStart)^2)+((dischargeXS1PointE-dischargeXS1EastStart)^2))
}

#To manually select ReferenceDistance:
dischargeXS1ReferenceDistance <- dischargePointsXS1$DistanceRaw[dischargePointsXS1$name=="DSC_LB_PIN"]


#Sets Horizontal adjustment value based on reference point coordinate.  
dischargeXS1HorizontalAdjust <- 0-dischargeXS1ReferenceDistance

#Transforms raw distance to adjusted distance based on reference distance point.
for(i in 1:(length(dischargePointsXS1$name))){
  dischargePointsXS1$DistanceAdj[i]<-dischargePointsXS1$DistanceRaw[i]+dischargeXS1HorizontalAdjust
}

# #Calculates the bankfull width
# DSCXS1Bankfull<-abs((dischargePointsXS1$DistanceAdj[grepl("RBF",dischargePointsXS1$name)])-
#                       (dischargePointsXS1$DistanceAdj[grepl("LBF",dischargePointsXS1$name)]))

#Creates dataframe of staff gauge points
staffGaugePoints=subset(surveyPtsDF,surveyPtsDF$name=="SP_0.55M")
staffGaugePoints<-staffGaugePoints[order(staffGaugePoints$N),]
rownames(staffGaugePoints)<-seq(length=nrow(staffGaugePoints))

#Set meter mark where the staff gauge was shot in and the name of the staff gauge point:
#Recorded in field data
staffGaugeMeterMark<-0.55
staffGaugeElevation <- staffGaugePoints$H[grepl("SP_0.55M",staffGaugePoints$name)]  

#Converts discharge XS1 transect point elevations to gauge height (rounded to 2 digits).
dischargePointsXS1$gaugeHeight<-dischargePointsXS1$H - (staffGaugeElevation - staffGaugeMeterMark)
dischargePointsXS1$gaugeHeight<-round(dischargePointsXS1$gaugeHeight,digits=2)

#Assigns a unique to each measurement for plot viewing purposes.  
dischargePointsXS1$ID<-c(1:length(dischargePointsXS1$name))

dischargePointsXS1 <- dischargePointsXS1[order(dischargePointsXS1$DistanceAdj),]

#Sets plot2 settings.  
xAxisTitle2<-list(title="Distance (m)",zeroline=FALSE, range=c(0,10))
yAxisTitle2<-list(title="Gauge Height  (m)",zeroline=FALSE)
font<-list(size=12,color='black')

#Plot the cross section by distance and gauge height.  
plot_ly(data=dischargePointsXS1,x=~DistanceAdj, y=~gaugeHeight, name='Distance vs. Gauge Height', type='scatter', mode='markers+lines', text=~name)%>%
  add_trace(y= 0,name = 'Gauge Height = 0.00m',mode='lines',line = list(color = 'red', width = 2, dash='dash')) %>%
  layout(title = siteID, xaxis=xAxisTitle2, yaxis=yAxisTitle2)

#####################################################################################################################################################
#Adjusts cross section elevations if necessary so lowest point is equal to 0.00m mark of staff gauge (prevents negative activation stage)
#####################################################################################################################################################
#Determines the lowest elevation of the discharge cross-section
dischargeXSmin<-min(dischargePointsXS1$H)

#Determines elevation of 0.00 meter mark of staff gage
staffGaugeZero=staffGaugeElevation-staffGaugeMeterMark

#Determines the offset between the lowest elevation and gauge height 
ElevOff<-dischargeXSmin-staffGaugeZero

#Adjusts the cross section elevations by the offset and rounds to 2 decimals
dischargePointsXS1$gaugeHeight<-dischargePointsXS1$gaugeHeight - ElevOff
dischargePointsXS1$gaugeHeight<-round(dischargePointsXS1$gaugeHeight,digits=2)

#Replots the adjusted cross section  
xAxisTitle2<-list(title="Distance (m)",zeroline=FALSE, range=c(0,10))
yAxisTitle2<-list(title="Gauge Height  (m)",zeroline=FALSE)
font<-list(size=12,color='black')
plot_ly(data=dischargePointsXS1,x=~DistanceAdj, y=~gaugeHeight, name='Distance vs. Gauge Height', type='scatter', mode='markers+lines', text=~name)%>%
  add_trace(y= 0,name = 'Gauge Height = 0.00m',mode='lines',line = list(color = 'red', width = 2, dash='dash')) %>%
  layout(title = siteID, xaxis=xAxisTitle2, yaxis=yAxisTitle2)
#####################################################################################################################################################

##### Now create the actual controls to upload... #####

#First, the addition or replacement when controls are activated table "geo_controlInfo_in"
numControls <- 2
geo_controlInfo_in_names <- c("locationID","startDate","endDate","controlNumber","segmentNumber","controlActivationState")
geo_controlInfo_in <- data.frame(matrix(nrow = numControls*numControls, ncol = length(geo_controlInfo_in_names)))
names(geo_controlInfo_in) <- geo_controlInfo_in_names

geo_controlInfo_in$locationID <- siteID
geo_controlInfo_in$startDate <- surveyActiveDate
geo_controlInfo_in$endDate <- surveyActiveDate
geo_controlInfo_in$controlNumber <- rep(1:numControls,numControls)
geo_controlInfo_in <- geo_controlInfo_in[order(geo_controlInfo_in$controlNumber),]
geo_controlInfo_in$segmentNumber <- rep(1:numControls,numControls)

#Known control activation states
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==geo_controlInfo_in$segmentNumber] <- 1
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber>geo_controlInfo_in$segmentNumber] <- 0

#Setting control activation states that are user defined.
#Is control #1 still active when control #2 is activated? 1 = Yes
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==1&geo_controlInfo_in$segmentNumber==2] <- 1

#Second, create entries for "geo_controlType_in" table for control parameters
geo_controlType_in_names <- c("locationID",
                               "startDate",
                               "endDate",
                               "controlNumber",
                               "hydraulicControlType",
                               "controlLeft",
                               "controlRight",
                               "rectangularWidth",
                               "rectangularWidthUnc",
                               "triangularAngle",
                               "triangularAngleUnc",
                               "parabolaWidth",
                               "parabolaWidthUnc",
                               "parabolaHeight",
                               "parabolaHeightUnc",
                               "orificeArea",
                               "orificeAreaUnc",
                               "channelSlope",
                               "channelSlopeUnc",
                               "manningCoefficient",
                               "manningCoefficientUnc",
                               "stricklerCoefficient",
                               "stricklerCoefficientUnc")
geo_controlType_in <- data.frame(matrix(nrow = numControls, ncol = length(geo_controlType_in_names)))
names(geo_controlType_in) <- geo_controlType_in_names

geo_controlType_in$locationID <- siteID
geo_controlType_in$startDate <- surveyActiveDate
geo_controlType_in$endDate <- surveyActiveDate
geo_controlType_in$controlNumber <- 1:numControls

#Entries for Control #1
geo_controlType_in$hydraulicControlType[1] <- "Rectangular Channel"
geo_controlType_in$controlLeft[1] <- dischargePointsXS1$DistanceAdj[dischargePointsXS1$name == "DSC_LFW"]
geo_controlType_in$controlRight[1] <- dischargePointsXS1$DistanceAdj[dischargePointsXS1$name == "DSC_XS33"]
geo_controlType_in$rectangularWidth[1] <- geo_controlType_in$controlRight[1]-geo_controlType_in$controlLeft[1]
geo_controlType_in$rectangularWidthUnc[1] <- 2.0 #Combined uncertainty associated with survey and where actual control begins (1.0 m default)

#Entries for Control #2
geo_controlType_in$hydraulicControlType[2] <- "Rectangular Channel"
geo_controlType_in$controlLeft[2] <- dischargePointsXS1$DistanceAdj[dischargePointsXS1$name == "DSC_XS33"]
geo_controlType_in$controlRight[2] <- dischargePointsXS1$DistanceAdj[dischargePointsXS1$name == "DSC_RBF"]
geo_controlType_in$rectangularWidth[2] <- geo_controlType_in$controlRight[2]-geo_controlType_in$controlLeft[2]
geo_controlType_in$rectangularWidthUnc[2] <- 2.0 #Combined uncertainty associated with survey and where actual control begins (1.0 m default)

#Slope calculations
colfunc <- colorRampPalette(c("cyan","deeppink"))
wettedEdgePoints=subset(surveyPtsDF,surveyPtsDF$mapCode%in%c("LEW","REW"))
wettedEdgePoints<-wettedEdgePoints[order(wettedEdgePoints$N),]
rownames(wettedEdgePoints)<-seq(length=nrow(wettedEdgePoints)) 
invisible(dev.new(noRStudioGD = TRUE))
# plot(wettedEdgePoints$E,wettedEdgePoints$N,pch=19, col=colfunc(length(wettedEdgePoints$H))[order(wettedEdgePoints$H)],
#      main=paste(siteID,"\nSelect a point above and below the discharge cross-section"),xlab="Raw Easting",ylab="Raw Northing")
# legend(min(wettedEdgePoints$E),max(wettedEdgePoints$N),legend=c("highest elevation","lowest elevation","discharge cross-section"),col = c("deeppink","cyan","green"),bty="n",pch = c(19,19,1))
# points(dischargePointsXS1$E,dischargePointsXS1$N, col="green")
# ans <- identify(wettedEdgePoints$E,wettedEdgePoints$N, n = 2, pos = F, tolerance = 0.25)
# Sys.sleep(1)
# invisible(dev.off())

ans=c(252,319)

#Plot subsetted wetted edges by manually entering ans values for tracking
wettedEdgePoints <- wettedEdgePoints[ans[1]:ans[2],]
# invisible(dev.new(noRStudioGD = TRUE))
# plot(wettedEdgePoints$E,wettedEdgePoints$N,pch=19, col=colfunc(length(wettedEdgePoints$H))[order(wettedEdgePoints$H)],
#      main=paste(siteID,"\nSelect two points above and below the discharge cross-section"),xlab="Raw Easting",ylab="Raw Northing")
# legend(min(wettedEdgePoints$E),max(wettedEdgePoints$N),legend=c("highest elevation","lowest elevation","discharge cross-section"),col = c("deeppink","cyan","green"),bty="n",pch = c(19,19,1))
# points(dischargePointsXS1$E,dischargePointsXS1$N, col="green")
# csOne <- identify(wettedEdgePoints$E,wettedEdgePoints$N, n = 2, pos = F, tolerance = 0.1)
# csTwo <- identify(wettedEdgePoints$E,wettedEdgePoints$N, n = 2, pos = F, tolerance = 0.1)
# Sys.sleep(1)
# invisible(dev.off())

csOne=c(9,7)
csTwo=c(54,56)

rise <- abs(mean(wettedEdgePoints$H[csOne])-mean(wettedEdgePoints$H[csTwo]))
run <- sqrt((mean(wettedEdgePoints$E[csOne])-mean(wettedEdgePoints$E[csTwo]))**2+(mean(wettedEdgePoints$N[csOne])-mean(wettedEdgePoints$N[csTwo]))**2)
geo_controlType_in$channelSlope[1] <- rise/run
geo_controlType_in$channelSlopeUnc[1] <- 0.005  #Default slope uncertainty is equal to slope
geo_controlType_in$channelSlope[2] <- rise/run
geo_controlType_in$channelSlopeUnc[2] <- 0.005 # Default slope uncertainty is equal to slope

#chosen to represent stream conditions with higher roughness above bankfull
geo_controlType_in$manningCoefficient[1] <- 0.05 # Cobble stream with some pools 
geo_controlType_in$manningCoefficientUnc[1] <- 0.025 # Default Mannings uncertainty equal 50%
geo_controlType_in$stricklerCoefficient[1] <- 1/geo_controlType_in$manningCoefficient[1]
geo_controlType_in$stricklerCoefficientUnc[1] <- geo_controlType_in$stricklerCoefficient[1]*(geo_controlType_in$manningCoefficientUnc[1]/geo_controlType_in$manningCoefficient[1])

geo_controlType_in$manningCoefficient[2] <- 0.1 # Trees and some brush 
geo_controlType_in$manningCoefficientUnc[2] <- 0.05 # Default Mannings uncertainty equal 50%
geo_controlType_in$stricklerCoefficient[2] <- 1/geo_controlType_in$manningCoefficient[2]
geo_controlType_in$stricklerCoefficientUnc[2] <- geo_controlType_in$stricklerCoefficient[2]*(geo_controlType_in$manningCoefficientUnc[2]/geo_controlType_in$manningCoefficient[2])

#Third,  use equations to populate "geo_priorParameters_in" table
geo_priorParameters_in <- data.frame(matrix(nrow = numControls, ncol = 10))
names(geo_priorParameters_in) <- c("locationID",
                                   "startDate",
                                   "endDate",
                                   "controlNumber",
                                   "priorExponent",
                                   "priorExponentUnc",
                                   "priorCoefficient",
                                   "priorCoefficientUnc",
                                   "priorActivationStage",
                                   "priorActivationStageUnc")

#Manually enter activation stages for controls
geo_priorParameters_in$priorActivationStage[1] <- dischargePointsXS1$gaugeHeight[dischargePointsXS1$name == "DSC_XS19"]
geo_priorParameters_in$priorActivationStageUnc[1] <- 0.2 # Combined uncertainty associated with survey and actual activation stage (0.1 m default)

geo_priorParameters_in$priorActivationStage[2] <- dischargePointsXS1$gaugeHeight[dischargePointsXS1$name == "DSC_XS43"]
geo_priorParameters_in$priorActivationStageUnc[2] <- 0.2 # Combined uncertainty associated with survey and actual activation stage (0.1 m default)

geo_priorParameters_in$locationID <- siteID
geo_priorParameters_in$startDate <- surveyActiveDate
geo_priorParameters_in$endDate <- surveyActiveDate

#Loop through to calculate exponent and coefficients
for(i in 1:numControls){
  geo_priorParameters_in$controlNumber[i] <- i
  if(!geo_controlType_in$hydraulicControlType[i]%in%c("Rectangular Weir","Rectangular Channel")){
    stop("Control type not found in the list.")
  }
  switch(geo_controlType_in$hydraulicControlType[i],
         "Rectangular Weir" = {
           Cr <- 0.4
           Cr_unc <- 0.1
           Bw <- geo_controlType_in$rectangularWidth[geo_controlType_in$controlNumber == i] #meters wide
           Bw_unc <- geo_controlType_in$rectangularWidthUnc[geo_controlType_in$controlNumber == i]
           g <- 9.81 #metersPerSecondSquared
           g_unc <- 0.01
           geo_priorParameters_in$priorCoefficient[i] <- Cr * Bw * (2*g)**(1/2)
           geo_priorParameters_in$priorCoefficientUnc[i] <- geo_priorParameters_in$priorCoefficient[i] * ((Cr_unc/Cr)**2+(Bw_unc/Bw)**2)**(1/2) + 0.5*(g_unc)/(2*g)
           
           geo_priorParameters_in$priorExponent[i] <- 1.5 #Recommended by BaM
           geo_priorParameters_in$priorExponentUnc[i] <- 0.05 #Recommended by BaM
         },
         "Rectangular Channel" = {
           Ks <- geo_controlType_in$stricklerCoefficient[geo_controlType_in$controlNumber == i]
           Ks_unc <- geo_controlType_in$stricklerCoefficientUnc[geo_controlType_in$controlNumber == i]
           Bw <- geo_controlType_in$rectangularWidth[geo_controlType_in$controlNumber == i] #meters wide
           Bw_unc <- geo_controlType_in$rectangularWidthUnc[geo_controlType_in$controlNumber == i]
           slope <- geo_controlType_in$channelSlope[geo_controlType_in$controlNumber == i]
           slope_unc <- geo_controlType_in$channelSlopeUnc[geo_controlType_in$controlNumber == i]
           geo_priorParameters_in$priorCoefficient[i] <- Ks * Bw * (slope)**(1/2)
           geo_priorParameters_in$priorCoefficientUnc[i] <- geo_priorParameters_in$priorCoefficient[i] * ((Ks_unc/Ks)**2+(Bw_unc/Bw)**2)**(1/2) + 0.5*(slope_unc)/(slope)
           
           geo_priorParameters_in$priorExponent[i] <- 1.67 #Recommended by BaM
           geo_priorParameters_in$priorExponentUnc[i] <- 0.05 #Recommended by BaM
         }
  )
}

#Plot controls to double check
invisible(dev.new(noRStudioGD = TRUE))
plot(dischargePointsXS1$DistanceAdj,dischargePointsXS1$gaugeHeight,main=paste(siteID,"Discharge XS1: Distance vs. Gauge Height"),xlab="Distance (m)",ylab="Gauge Height (m)")
text(dischargePointsXS1$DistanceAdj,dischargePointsXS1$gaugeHeight,labels=dischargePointsXS1$name,pos=4)
lines(lines(dischargePointsXS1$DistanceAdj,dischargePointsXS1$gaugeHeight,lty=3))
colorsForPlot <- c("blue","red","green","orange","purple")
for(i in 1:numControls){
  x <- c(geo_controlType_in$controlLeft[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlLeft[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlRight[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlRight[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlLeft[geo_controlType_in$controlNumber==i])
  
  #Determine ymax
  if(i == numControls){
    ymax <- max(dischargePointsXS1$gaugeHeight)
  }else if(any(geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==i&geo_controlInfo_in$segmentNumber>i]==0)){
    overtakingControlNumber <- min(geo_controlInfo_in$segmentNumber[geo_controlInfo_in$controlNumber==i&
                                                                      geo_controlInfo_in$segmentNumber>i&
                                                                      geo_controlInfo_in$controlActivationState==0])
    ymax <- geo_priorParameters_in$priorActivationStage[geo_priorParameters_in$controlNumber == overtakingControlNumber]
  }else{
    ymax <- max(dischargePointsXS1$gaugeHeight)
  }
  
  #Determine ymin if control overtakes others
  if(i == 1){
    ymin <- geo_priorParameters_in$priorActivationStage[geo_priorParameters_in$controlNumber==i]
  }else if(sum(geo_controlInfo_in$controlActivationState[geo_controlInfo_in$segmentNumber==i&geo_controlInfo_in$controlNumber<i])<(i-1)){
    #ymin <- min(geo_priorParameters_in$priorActivationStage[geo_priorParameters_in$controlNumber<i])
    ymin <- geo_priorParameters_in$priorActivationStage[geo_priorParameters_in$controlNumber==i]
  }else{
    ymin <- geo_priorParameters_in$priorActivationStage[geo_priorParameters_in$controlNumber==i]
  }
  
  y <- c(ymax,
         ymin,
         ymin,
         ymax,
         ymax)
  polygon(x,y, col = adjustcolor(colorsForPlot[i],alpha.f = 0.5))
}
dev.copy2pdf(file = paste0(filepath,"/",siteID,"_siteControls.pdf"), width = 16, height = 9)

#Write out three tables for ingest to GitHub and testing location both
geo_controlInfo_in_output <- c("locationID",
                               "startDate",
                               "endDate",
                               "controlNumber",
                               "segmentNumber",
                               "controlActivationState")
geo_controlInfo_in <- geo_controlInfo_in[,names(geo_controlInfo_in)%in%geo_controlInfo_in_output]
write.csv(geo_controlInfo_in,
          paste0(filepath,"/geo_controlInfo_in.csv"),
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")

geo_controlType_in_output <- c("locationID",
                               "startDate",
                               "endDate",
                               "controlNumber",
                               "hydraulicControlType",
                               "rectangularWidth",
                               "rectangularWidthUnc",
                               "triangularAngle",
                               "triangularAngleUnc",
                               "parabolaWidth",
                               "parabolaWidthUnc",
                               "parabolaHeight",
                               "parabolaHeightUnc",
                               "orificeArea",
                               "orificeAreaUnc",
                               "channelSlope",
                               "channelSlopeUnc",
                               "manningCoefficient",
                               "manningCoefficientUnc",
                               "stricklerCoefficient",
                               "stricklerCoefficientUnc")
geo_controlType_in <- geo_controlType_in[,names(geo_controlType_in)%in%geo_controlType_in_output]
write.csv(geo_controlType_in,
          paste0(filepath,"/geo_controlType_in.csv"),
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")

geo_priorParameters_in_output <- c("locationID",
                                   "startDate",
                                   "endDate",
                                   "controlNumber",
                                   "priorExponent",
                                   "priorExponentUnc",
                                   "priorCoefficient",
                                   "priorCoefficientUnc",
                                   "priorActivationStage",
                                   "priorActivationStageUnc")
geo_priorParameters_in <- geo_priorParameters_in[,names(geo_priorParameters_in)%in%geo_priorParameters_in_output]
write.csv(geo_priorParameters_in,
          paste0(filepath,"/geo_priorParameters_in.csv"),
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")

