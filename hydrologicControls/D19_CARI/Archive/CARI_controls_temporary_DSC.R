
# library(neonUtilities)
# #For use with the API functionality
# dataDir <- "API"
# siteID <- "HOPB"
# 
# streamMorphoDPID <- "DP4.00131.001"
# dataFromAPI <- zipsByProduct(streamMorphoDPID,site,package="expanded",check.size=FALSE)
# filepath <- "C:/Users/kcawley/Desktop/test/filesToStack00131/"
# stackByTable(filepath=filepath, folder = TRUE)

#Use the new function here once we have the data in the new ECS zip packages

#HOPB processing code
# filepath <- "C:/Users/kcawley/Desktop/test/"
# surveyPtsDF <- read.table(paste0(filepath,"HOPB_surveyPts_20170921.csv"),sep = ",",stringsAsFactors = FALSE, header = TRUE)
# siteID <- "HOPB"
# surveyDate <- "2017-09-21T16:00"

#load packages
require(rgdal)
library(rgdal)
library(plotly)

#NEON Domain number (ex: D01).
domainID<-'D19' 

#Four-digit NEON site code (ex: HOPB).
siteID <- 'CARI'  

#The end date of the geomorphology survey (YYYYMMDD).  
surveyDate<-'20170913' 

#Stipulate 4-digit site code, underscore, and survey year (ex: HOPB_2017). 
surveyID <- "CARI_2017"  

#Queues a directory that contains file paths for each site per survey date.  
siteDirectory<-read.csv('N:/Science/AQU/Geomorphology_Survey_Data/inputDirectory.csv',head=T,sep=",",stringsAsFactors = F) 

#Creates dataframe that contains survey data.  
filePath <- siteDirectory$filePath[which(siteDirectory$surveyID==surveyID)]
surveyShapefileName <- siteDirectory$surveyShapefileName[which(siteDirectory$surveyID==surveyID)]
surveyPts <- readOGR(filePath,surveyShapefileName)
surveyPtsDF <- as.data.frame(surveyPts)

#Working directory where files will be output.  
wdir<-paste('C:/Users/nharrison/Documents/GitHub/landWaterSoilIPT/streamMorpho/ScienceProcessingCode/R_Metrics',siteID,'Raw_Data',sep="/") 
#wdir<-paste('C:/Users/kcawley/Documents/GitHub/landWaterSoilIPT/streamMorpho/ScienceProcessingCode/R_Metrics',siteID,'Raw_Data',sep="/") 

#Creates dataframe of all points associated with transect DSC2 (the temporary discharge transect).  
dischargePointsXS2<-subset(surveyPtsDF,mapCode=="Transect_DSC2")
dischargePointsXS2<-dischargePointsXS2[order(-dischargePointsXS2$E),]
rownames(dischargePointsXS2)<-seq(length=nrow(dischargePointsXS2))

#Sets plot1 settings.  
xAxisTitle1<-list(title="Easting (m)",zeroline=FALSE)
yAxisTitle1<-list(title="Northing  (m)",zeroline=FALSE)
font<-list(size=12,color='black')

#Plot the cross section by easting and northing data.
plot_ly(data=dischargePointsXS2,x=~E, y=~N, name='Easting vs Northing', type='scatter', mode='markers', text=~name)%>%
  layout(title = siteID, xaxis=xAxisTitle1, yaxis=yAxisTitle1)

#Manually select NorthStart and EastStart coordinates
dischargeXS2NorthStart<-dischargePointsXS2$N[25]
dischargeXS2EastStart<-dischargePointsXS2$E[25]

#Assigns a raw Distance value to each point relative to the NorthStart and EastStart coordinates.
for(i in 1:(length(dischargePointsXS2$name))){
  dischargeXS2PointN<-dischargePointsXS2$N[i]
  dischargeXS2PointE<-dischargePointsXS2$E[i]
  dischargePointsXS2$DistanceRaw[i]<-sqrt(((dischargeXS2PointN-dischargeXS2NorthStart)^2)+((dischargeXS2PointE-dischargeXS2EastStart)^2))
}

#To manually select ReferenceDistance:
dischargeXS2ReferenceDistance<-dischargePointsXS2$DistanceRaw[25]

#Sets Horizontal adjustment value based on reference point coordinate.  
dischargeXS2HorizontalAdjust<-0-dischargeXS2ReferenceDistance

#Transforms raw distance to adjusted distance based on reference distance point.
for(i in 1:(length(dischargePointsXS2$name))){
  dischargePointsXS2$DistanceAdj[i]<-dischargePointsXS2$DistanceRaw[i]+dischargeXS2HorizontalAdjust
}

#Plot the cross-section again to check reference values and ensure that the profile is being viewed from left to right bank.

#Sets plot2 settings.  
xAxisTitle2<-list(title="Distance (m)",zeroline=FALSE)
yAxisTitle2<-list(title="Elevation  (m)",zeroline=FALSE)
font<-list(size=12,color='black')

#Plot the cross section by easting and northing data.
plot_ly(data=dischargePointsXS2,x=~DistanceAdj, y=~H, name='Distance vs Elevation', type='scatter', mode='lines+markers', text=~name)%>%
  layout(title = siteID, xaxis=xAxisTitle2, yaxis=yAxisTitle2)

#Calculates the bankfull width.
DSCXS2Bankfull<-abs((dischargePointsXS2$DistanceAdj[grepl("RBF",dischargePointsXS2$name)])-(dischargePointsXS2$DistanceAdj[grepl("LBF",dischargePointsXS2$name)]))

#Creates dataframe of staff gauge points.
staffGaugePoints=subset(surveyPtsDF,surveyPtsDF$mapCode=="Gauge_Transect_DSC2")
staffGaugePoints<-staffGaugePoints[order(staffGaugePoints$N),]
rownames(staffGaugePoints)<-seq(length=nrow(staffGaugePoints))

#Set meter mark where the staff gauge was shot in and the name of the staff gauge point:
#Recorded in field data
staffGaugeMeterMark<-0.40
staffGaugeElevation <- staffGaugePoints$H[grepl("SP_0.40M_2",staffGaugePoints$name)]  

#Converts discharge XS1 transect point elevations to gauge height (rounded to 2 digits).
dischargePointsXS2$gaugeHeight<-dischargePointsXS2$H - (staffGaugeElevation - staffGaugeMeterMark)
dischargePointsXS2$gaugeHeight<-round(dischargePointsXS2$gaugeHeight,digits=2)

#Assigns a unique to each measurement for plot viewing purposes.  
dischargePointsXS2$ID<-c(1:length(dischargePointsXS2$name))

dischargePointsXS2 <- dischargePointsXS2[order(dischargePointsXS2$DistanceAdj),]
#invisible(dev.new(noRStudioGD = TRUE))

#Sets plot3 settings.  
xAxisTitle3<-list(title="Distance (m)",zeroline=FALSE) #Define range , range=c(-5,15)
yAxisTitle3<-list(title="Gauge Height  (m)",zeroline=FALSE)
font<-list(size=12,color='black')

#Plot the cross section by distance and gauge height.  Note whether or not red line is below thalweg.  
plot_ly(data=dischargePointsXS2,x=~DistanceAdj, y=~gaugeHeight, name='Distance vs. Gauge Height', type='scatter', mode='markers+lines', text=~name)%>%
  add_trace(y= 0,name = 'Gauge Height = 0.00m',mode='lines',line = list(color = 'red', width = 2, dash='dash')) %>%
  layout(title = siteID, xaxis=xAxisTitle3, yaxis=yAxisTitle3)

#Calculates gaugeHeight at LB and RB bankfull:
gaugeHeightLBF<-dischargePointsXS2$gaugeHeight[grepl("LBF",dischargePointsXS2$name)]
gaugeHeightRBF<-dischargePointsXS2$gaugeHeight[grepl("RBF",dischargePointsXS2$name)]

#Asseses whether negative stage is present in the discharge cross-section 
negativeStage<-any(dischargePointsXS2$gaugeHeight<0)

if(negativeStage==TRUE){
  exectpart<-TRUE
}else{exectpart<-FALSE}

#This section will only run if there are negative values in the gauge height column.  
if(exectpart){
  
  #Determines the lowest elevation of the discharge cross-section, assumed to be the thalweg. 
  dischargeXS2THL<-min(dischargePointsXS2$H)
  
  #Calculates the elevation of the 0.00 meter mark of the staff gauge.  
  staffGaugeZeroElevation<-(staffGaugePoints$H[staffGaugePoints$name=="SP_0.40M_2"])-staffGaugeMeterMark
  
  #Calculates the difference between the staff gauge 0.00m mark elevation and the discharge thalweg elevation.   
  gaugeZeroQElevDiff<--as.numeric(dischargeXS2THL-staffGaugeZeroElevation)
  
  #Offsets the elevation of the gauge heights by this difference, rounds to two digits.  
  dischargePointsXS2$gaugeOffsetElevation<-dischargePointsXS2$H + (gaugeZeroQElevDiff)
  dischargePointsXS2$gaugeOffsetElevation<-round(dischargePointsXS2$gaugeOffsetElevation,digits=2)
  
  #Offsets the gauge heights by the offset elevation, rounds to two digits.  
  dischargePointsXS2$gaugeHeightOffset<-dischargePointsXS2$gaugeOffsetElevation - (staffGaugeElevation - staffGaugeMeterMark)
  dischargePointsXS2$gaugeHeightOffset<-round(dischargePointsXS2$gaugeHeightOffset,digits=2)
  
  #Plots discharge XS1 transect point distances vs gaugeHeightOffset.  Red line should be at the thalweg.
  plot_ly(data=dischargePointsXS2,x=~DistanceAdj, y=~gaugeHeightOffset, name='Distance vs. Gauge Height', type='scatter', mode='markers+lines', text=~name)%>%
    add_trace(y= 0,name = 'Gauge Height = 0.00m',mode='lines',line = list(color = 'red', width = 2, dash='dash')) %>%
    layout(title = paste(siteID,":","Temporary DSC XS, Distance vs. Gauge Height"), xaxis=xAxisTitle3, yaxis=yAxisTitle3)
  
  
}else{"There are no negative gauge height values in discharge XS2.  There is no need for correction."}

#### Now create the actual controls to upload... #####

#First, the addition or replacement when controls are activated table "geo_controlInfo_in"
numControls <- 3
geo_controlInfo_in_names <- c("locationID","startDate","endDate","controlNumber","segmentNumber","controlActivationState")
geo_controlInfo_in <- data.frame(matrix(nrow = numControls*numControls, ncol = length(geo_controlInfo_in_names)))
names(geo_controlInfo_in) <- geo_controlInfo_in_names

geo_controlInfo_in$locationID <- siteID
geo_controlInfo_in$startDate <- surveyDate
geo_controlInfo_in$endDate <- surveyDate
geo_controlInfo_in$controlNumber <- rep(1:numControls,numControls)
geo_controlInfo_in <- geo_controlInfo_in[order(geo_controlInfo_in$controlNumber),]
geo_controlInfo_in$segmentNumber <- rep(1:numControls,numControls)

#Known control activation states
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==geo_controlInfo_in$segmentNumber] <- 1
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber>geo_controlInfo_in$segmentNumber] <- 0

#Setting control activation states that are user defined.

#Is control #1 still active when control #2 is activated? 1 = Yes
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==1&geo_controlInfo_in$segmentNumber==2] <- 1

#Is control #1 still active when control #3 is activated? 0 = No
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==1&geo_controlInfo_in$segmentNumber==3] <- 0

#Is control #2 still active when control #3 is activated? 0 = No
geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==2&geo_controlInfo_in$segmentNumber==3] <- 0

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
geo_controlType_in$startDate <- surveyDate
geo_controlType_in$endDate <- surveyDate
geo_controlType_in$controlNumber <- 1:numControls

#Entries for Control #1
geo_controlType_in$hydraulicControlType[1] <- "Rectangular Weir"
geo_controlType_in$controlLeft[1] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "18"]
# geo_controlType_in$controlLeft[1] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$name == "DSC21"]
geo_controlType_in$controlRight[1] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "41"]
# geo_controlType_in$controlRight[1] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$name == "DSC29"]
geo_controlType_in$rectangularWidth[1] <- geo_controlType_in$controlRight[1]-geo_controlType_in$controlLeft[1]
geo_controlType_in$rectangularWidthUnc[1] <- 0.05 #Uncertainty associated with AIS survey

#Entries for Control #2
geo_controlType_in$hydraulicControlType[2] <- "Rectangular Weir"
geo_controlType_in$controlLeft[2] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "13"]
geo_controlType_in$controlRight[2] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "18"]
# geo_controlType_in$controlLeft[2] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$name == "DSC8"]
# geo_controlType_in$controlRight[2] <- dischargePointsXS2$DistanceAdj[dischargePointsXS2$name == "DSC34"]
geo_controlType_in$rectangularWidth[2] <- geo_controlType_in$controlRight[2]-geo_controlType_in$controlLeft[2]
geo_controlType_in$rectangularWidthUnc[2] <- 0.05

#Entries for Control #3
geo_controlType_in$hydraulicControlType[3] <- "Rectangular Channel"
geo_controlType_in$controlLeft[3] <- (dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "1"]+dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "4"])/2
geo_controlType_in$controlRight[3] <- (dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "41"]+dischargePointsXS2$DistanceAdj[dischargePointsXS2$ID == "43"])/2
geo_controlType_in$rectangularWidth[3] <- geo_controlType_in$controlRight[3]-geo_controlType_in$controlLeft[3]
geo_controlType_in$rectangularWidthUnc[3] <- 0.05

#Slope calculations
colfunc <- colorRampPalette(c("cyan","deeppink"))
wettedEdgePoints=subset(surveyPtsDF,surveyPtsDF$mapCode%in%c("LEW","REW"))
wettedEdgePoints<-wettedEdgePoints[order(wettedEdgePoints$N),]
rownames(wettedEdgePoints)<-seq(length=nrow(wettedEdgePoints)) 
invisible(dev.new(noRStudioGD = TRUE))
plot(wettedEdgePoints$E,wettedEdgePoints$N,pch=19, col=colfunc(length(wettedEdgePoints$H))[order(wettedEdgePoints$H)],
     main=paste(siteID,"\nSelect a point above and below the discharge cross-section"),xlab="Raw Easting",ylab="Raw Northing")
legend(min(wettedEdgePoints$E),max(wettedEdgePoints$N),legend=c("highest elevation","lowest elevation","discharge cross-section"),col = c("deeppink","cyan","green"),bty="n",pch = c(19,19,1))
points(dischargePointsXS2$E,dischargePointsXS2$N, col="green")
ans <- identify(wettedEdgePoints$E,wettedEdgePoints$N, n = 2, pos = F, tolerance = 0.25)
#ans = 981, 932
Sys.sleep(1)
invisible(dev.off())

#Plot subsetted wetted edges by manually entering ans values for tracking
wettedEdgePoints <- wettedEdgePoints[932:981,]
invisible(dev.new(noRStudioGD = TRUE))
plot(wettedEdgePoints$E,wettedEdgePoints$N,pch=19, col=colfunc(length(wettedEdgePoints$H))[order(wettedEdgePoints$H)],
     main=paste(siteID,"\nSelect two points above and below the discharge cross-section"),xlab="Raw Easting",ylab="Raw Northing")
legend(min(wettedEdgePoints$E),max(wettedEdgePoints$N),legend=c("highest elevation","lowest elevation","discharge cross-section"),col = c("deeppink","cyan","green"),bty="n",pch = c(19,19,1))
points(dischargePointsXS2$E,dischargePointsXS2$N, col="green")
csOne <- identify(wettedEdgePoints$E,wettedEdgePoints$N, n = 2, pos = F, tolerance = 0.1)
#csOne = 22,23
csTwo <- identify(wettedEdgePoints$E,wettedEdgePoints$N, n = 2, pos = F, tolerance = 0.1)
#csTwo = 14,15
Sys.sleep(1)
invisible(dev.off())

rise <- abs(mean(wettedEdgePoints$H[csOne])-mean(wettedEdgePoints$H[csTwo]))
run <- sqrt((mean(wettedEdgePoints$E[csOne])-mean(wettedEdgePoints$E[csTwo]))**2+(mean(wettedEdgePoints$N[csOne])-mean(wettedEdgePoints$N[csTwo]))**2)
geo_controlType_in$channelSlope[3] <- rise/run
geo_controlType_in$channelSlopeUnc[3] <- 0.015

#chosen to represent stream conditions with higher roughness above bankfull
geo_controlType_in$manningCoefficient[3] <- 0.05
geo_controlType_in$manningCoefficientUnc[3] <- 0.001
geo_controlType_in$stricklerCoefficient[3] <- 1/geo_controlType_in$manningCoefficient[3]
geo_controlType_in$stricklerCoefficientUnc[3] <- geo_controlType_in$stricklerCoefficient[3]*(geo_controlType_in$manningCoefficientUnc[3]/geo_controlType_in$manningCoefficient[3])

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
geo_priorParameters_in$priorActivationStage[1] <- dischargePointsXS2$gaugeHeight[dischargePointsXS2$ID == "29"]
geo_priorParameters_in$priorActivationStageUnc[1] <- 0.01

geo_priorParameters_in$priorActivationStage[2] <- dischargePointsXS2$gaugeHeight[dischargePointsXS2$ID == "18"]
geo_priorParameters_in$priorActivationStageUnc[2] <- 0.01

geo_priorParameters_in$priorActivationStage[3] <- (dischargePointsXS2$gaugeHeight[dischargePointsXS2$ID == "13"]+dischargePointsXS2$gaugeHeight[dischargePointsXS2$ID == "41"])/2
geo_priorParameters_in$priorActivationStageUnc[3] <- 0.01

geo_priorParameters_in$locationID <- siteID
geo_priorParameters_in$startDate <- surveyDate
geo_priorParameters_in$endDate <- surveyDate

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
plot(dischargePointsXS2$DistanceAdj,dischargePointsXS2$gaugeHeight,main=paste(siteID,"Discharge XS1: Distance vs. Gauge Height"),xlab="Distance (m)",ylab="Gauge Height (m)")
text(dischargePointsXS2$DistanceAdj,dischargePointsXS2$gaugeHeight,labels=dischargePointsXS2$name,pos=4)
lines(lines(dischargePointsXS2$DistanceAdj,dischargePointsXS2$gaugeHeight,lty=3))
colorsForPlot <- c("blue","red","green","orange","purple")
for(i in 1:numControls){
  x <- c(geo_controlType_in$controlLeft[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlLeft[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlRight[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlRight[geo_controlType_in$controlNumber==i],
         geo_controlType_in$controlLeft[geo_controlType_in$controlNumber==i])
  
  #Determine ymax
  if(i == numControls){
    ymax <- max(dischargePointsXS2$gaugeHeight)
  }else if(any(geo_controlInfo_in$controlActivationState[geo_controlInfo_in$controlNumber==i&geo_controlInfo_in$segmentNumber>i]==0)){
    overtakingControlNumber <- min(geo_controlInfo_in$segmentNumber[geo_controlInfo_in$controlNumber==i&
                                                                      geo_controlInfo_in$segmentNumber>i&
                                                                      geo_controlInfo_in$controlActivationState==0])
    ymax <- geo_priorParameters_in$priorActivationStage[geo_priorParameters_in$controlNumber == overtakingControlNumber]
  }else{
    ymax <- max(dischargePointsXS2$gaugeHeight)
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
dev.copy2pdf(file = paste0(siteID,"_siteControls.pdf"), width = 16, height = 9)

#Write out three tables for ingest to GitHub and testing location both
geo_controlInfo_in_output <- c("locationID",
                               "startDate",
                               "endDate",
                               "controlNumber",
                               "segmentNumber",
                               "controlActivationState")
geo_controlInfo_in <- geo_controlInfo_in[,names(geo_controlInfo_in)%in%geo_controlInfo_in_output]
write.csv(geo_controlInfo_in,
          "geo_controlInfo_in.csv",
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(geo_controlInfo_in,
          "H:/controlTesting/geo_controlInfo_in.csv",
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
          "geo_controlType_in.csv",
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(geo_controlType_in,
          "H:/controlTesting/geo_controlType_in.csv",
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
          "geo_priorParameters_in.csv",
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")
write.csv(geo_priorParameters_in,
          "H:/controlTesting/geo_priorParameters_in.csv",
          quote = TRUE,
          row.names = FALSE,
          fileEncoding = "UTF-8")
