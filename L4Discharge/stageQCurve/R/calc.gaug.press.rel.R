##############################################################################################
#' @title Function to calculate Gauge-Pressure relationship

#' @author
#' Zachary Nickerson \email{nickerson@battelleecology.org} \cr

#' @description This script uses multiple data products to establish a relationship between
#' stream stage and calibrated pressure for a given water year

#' @param DIRPATH An environment variable that contains the location of the files in
#' the Docker container [string]
#' @param BAMWS An environment variable that contains the location of the BaM config
#' files in the Docker container [string]
#' @param DATAWS An environment variable that contains the location of the downloaded NEON
#' zip files in the Docker container [string]
#' @param startDate A date in the format of "YYYY-MM-DD", entered as env variable in
#' run command [string]
#' @param site A four letter site code for a NEON site, entered as env variable in
#' run command [string]

#' @return This function results in a text file being written to DATAWS.

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @export

# changelog and author contributions / copyrights
#   Zachary Nickerson (2023-07-25)
#     original creation
##############################################################################################
calc.gaug.press.rel <- function(DIRPATH = Sys.getenv("DIRPATH"),
                                BAMWS = Sys.getenv("BAMWS"),
                                DATAWS = Sys.getenv("DATAWS"),
                                startDate = Sys.getenv("STARTDATE"),
                                site = Sys.getenv("SITE")){
  
  # String constants
  osDateFormat <- "%Y-%m-%dT%H:%M:%S"
  osPubDateFormat <- "%Y-%m-%dT%H:%MZ"
  isDateFormat <- "%Y-%m-%dT%H:%M:%S.000Z"
  
  # Numeric Constants
  secondsInMinute <- 60
  secondsInDay <- 60*60*24
  secondsIn10Min <- 10*60
  secondsInYear <- 60*60*24*365
  minutesInDay <- 60*24
  #gaugeSurveyLoc <- 0.8#meters
  convKPatoPa <- 1000 #Pa per 1 kPa
  roe = 999 #kg/m^3 density of water
  grav = 9.80665 #m/s^2 gravity constant
  
  # Now for when I need it
  now <- as.POSIXct(Sys.time())
  
  # Metadata
  searchIntervalStartDate <- as.POSIXct(startDate, tz = "UTC")
  WYsearchIntervalStartDate <- as.POSIXct(format(stageQCurve::def.calc.WY.strt.end.date(searchIntervalStartDate)[["startDate"]], format = "%Y-%m-%d"))
  WYsearchIntervalEndDate <- as.POSIXct(format(stageQCurve::def.calc.WY.strt.end.date(searchIntervalStartDate)[["endDate"]], format = "%Y-%m-%d"))+secondsInDay
  yearMonth <- format(searchIntervalStartDate,format="%Y-%m")
  year <- format(searchIntervalStartDate,format="%Y")
  month <- format(searchIntervalStartDate,format="%m")
  downloadedDataPath <- DATAWS
  if (Sys.getenv("SITE")%in%c("BLWA","TOMB","FLNT","TOOK")) {
    mod <- "bat"
  }else{
    mod <- "geo"
  }
  
  # Get continuous discharge data from the OS system
  # If data has been downloaded using neonUtilities::zipsByProduct() and saved to DATAWS
  if(file.exists(paste0(downloadedDataPath,"filesToStack00130"))){
    continuousData  <- try(read.csv(paste(downloadedDataPath,"filesToStack00130","stackedFiles","csd_continuousDischarge.csv", sep = "/")),silent = T)
    regressionData <- try(read.csv(paste(downloadedDataPath,"filesToStack00130","stackedFiles",paste0(mod,"_gaugeWaterColumnRegression.csv"), sep = "/")),silent = T)
  }else{
    stop("Unzipped and stacked data from DP4.00130.001 could not be retrieved. Follow instructions in the README of the package to properly format input data")
  }
  
  # Error handling if no discharge or curve identification data can be found
  if(attr(continuousData, "class") == "try-error"){
    failureMessage <- paste0("Data could not be retrieved from DP4.00130.001:csd_continuousDischarge_pub. Ensure the required input data are stored in ",downloadedDataPath)
    stop(failureMessage)
  }
  if(attr(regressionData, "class") == "try-error"){
    failureMessage <- paste0("Data could not be retrieved from DP4.00130.001:geo/bat_gaugeWaterColumnRegression. Ensure the required input data are stored in ",downloadedDataPath)
    stop(failureMessage)
  }
  
  # Subset to the metadata and calibrated pressure data for the time range specified
  continuousData$CalPress <- continuousData$calibratedPressure
  continuousData$validCalFlag <- continuousData$dischargeValidCalQF
  continuousData$date <- as.POSIXct(continuousData$endDate, format = osPubDateFormat)
  L0PressureData <- continuousData[names(continuousData)%in%c("siteID","namedLocation","stationHorizontalID","startDate","date","nonSystematicUnc","CalPress","validCalFlag","sWatElevFinalQFSciRvw")]
  
  # Error handling is there is no data
  if(nrow(L0PressureData)<1){
    failureMessage <- paste0("Zero (0) records available from DP4.00130.001 for ",site," ",year,"-",month)
    stop(failureMessage)
  }
  
  #Get location information for the gauge location at the site
  gaugeLocData <- try(geoNEON::getLocBySite(site,
                                            type = "AQU",
                                            history = T),
                      silent = T)
  # Error handling if the GET was unsuccessful
  if(attr(gaugeLocData, "class") == "try-error"){
    failureMessage <- "Gauge location history could not be retrieved"
    stop(failureMessage)
  }
  gaugeLocData <- gaugeLocData[grepl("gauge",gaugeLocData$namedLocation),]
  #Error handling if there is no data
  if(nrow(gaugeLocData)<1){
    failureMessage <- "Zero (0) gauge location history records were retrieved"
    stop(gaugeLocData)
  }
  #Set data types
  #Add an endDate to avoid NA and null situations
  if(is.null(gaugeLocData$locationEndDate)){
    gaugeLocData$locationEndDate[is.null(gaugeLocData$locationEndDate)] <- format(as.POSIXct(now),format = osDateFormat)}
  if(any(is.na(gaugeLocData$locationEndDate))){
    gaugeLocData$locationEndDate[is.na(gaugeLocData$locationEndDate)] <- format(as.POSIXct(now),format = osDateFormat)}
  gaugeLocData$locationStartDate <-  as.POSIXct(gaugeLocData$locationStartDate,format = osDateFormat)
  gaugeLocData$locationEndDate <- as.POSIXct(gaugeLocData$locationEndDate,format = osDateFormat)
  gaugeLocData$zOffset <- as.numeric(gaugeLocData$zOffset)
  gaugeLocData$elevation <- as.numeric(gaugeLocData$elevation)
  gaugeLocData <- gaugeLocData[order(gaugeLocData$locationEndDate),]
  gaugeLocData$refElevPlusZ <- gaugeLocData$elevation + gaugeLocData$zOffset
  gaugeNamedLocation <- unique(gaugeLocData$namedLocation)
  
  # Beginning here, TOOK transitions will loop through each location
  if(site=="TOOK"){
    numLocs <- 2
  }else{
    numLocs <- 1
  }
  
  for(l in 1:numLocs){
    if(site=="TOOK"){
      searchGaugeNamedLocation <- gaugeNamedLocation[grepl(c("inflow","outflow")[l],gaugeNamedLocation)]
      dischargeNamedLocation <- paste0(site,c(".AOS.discharge.inflow",".AOS.discharge.outflow")[l])
    }else{
      searchGaugeNamedLocation <- gaugeNamedLocation
      dischargeNamedLocation <- paste0(site,".AOS.discharge")
    }
    
    # Get all relevant gauge height information to develop relationship
    # Gauge height data can be published in 3 locations: DP1.20267.001:gag_fieldData, DP1.20048.001:dsc_fieldData, DP1.20048.001:dsc_fieldDataADCP
    # If data has been downloaded using neonUtilities::zipsByProduct() and saved to DATAWS
    if(file.exists(paste0(downloadedDataPath,"filesToStack20267"))){
      L1GaugeData  <- try(read.csv(paste(downloadedDataPath,"filesToStack20267","stackedFiles","gag_fieldData.csv", sep = "/")),silent = T)
    }else{
      stop("Unzipped and stacked data from DP1.20267.001 could not be retrieved. Follow instructions in the README of the package to properly format input data")
    }
    if(file.exists(paste0(downloadedDataPath,"filesToStack20048"))){
      # Discharge data may not be available from one of the two tables
      flowmeterDischargeData  <- suppressWarnings(try(read.csv(paste(downloadedDataPath,"filesToStack20048","stackedFiles","dsc_fieldData.csv", sep = "/")),silent = T))
      adcpDischargeData  <- suppressWarnings(try(read.csv(paste(downloadedDataPath,"filesToStack20048","stackedFiles","dsc_fieldDataADCP.csv", sep = "/")),silent = T))
    }else{
      stop("Unzipped and stacked data from DP1.20048.001 could not be retrieved. Follow instructions in the README of the package to properly format input data")
    }
    
    # Error handling if no gauge data can be found
    if(attr(L1GaugeData, "class") == "try-error"){
      failureMessage <- paste0("Data could not be retrieved from DP1.20267.001:gag_fieldData. Ensure the required input data are stored in ",downloadedDataPath)
      stop(failureMessage)
    }
    if((attr(flowmeterDischargeData, "class") == "try-error")&(attr(adcpDischargeData, "class") == "try-error")){
      failureMessage <- paste0("Data could not be retrieved from DP1.20048.001:dsc_fieldData and DP1.20048.001:dsc_fieldDataADCP. Ensure the required input data are stored in ",downloadedDataPath)
      stop(failureMessage)
    }else{
      # Format event IDs
      if(!(attr(flowmeterDischargeData, "class") == "try-error")){
        if (site=="TOOK") {
          if(grepl("inflow",dischargeNamedLocation)){
            flowmeterDischargeData$eventID <- paste("TKIN",format(as.Date(flowmeterDischargeData$collectDate),"%Y%m%d"),sep = ".")
          }else{
            if(grepl("outflow",dischargeNamedLocation)){
              flowmeterDischargeData$eventID <- paste("TKOT",format(as.Date(flowmeterDischargeData$collectDate),"%Y%m%d"),sep = ".")
            }
          }
        }else{
          flowmeterDischargeData$eventID <- paste(flowmeterDischargeData$siteID,format(as.Date(flowmeterDischargeData$collectDate),"%Y%m%d"),sep = ".")
        }
      }
      if(!(attr(adcpDischargeData, "class") == "try-error")){
        adcpDischargeData$eventID <- paste(adcpDischargeData$siteID,format(as.Date(adcpDischargeData$endDate),"%Y%m%d"),sep = ".")
        adcpDischargeData$finalDischarge <- as.numeric(adcpDischargeData$totalDischarge)*1000
      }
    }
    
    # Join all gauge data into a single data frame and assign gauge location data
    if(!(attr(flowmeterDischargeData, "class") == "try-error")&(attr(adcpDischargeData, "class") == "try-error")){
      stageData <- flowmeterDischargeData
    }else{
      if((attr(flowmeterDischargeData, "class") == "try-error")&!(attr(adcpDischargeData, "class") == "try-error")){
        stageData <- adcpDischargeData
      }else{
        stageData <- merge(flowmeterDischargeData,adcpDischargeData,all = T)
      }
    }
    stageData$initialStageHeight <- stageData$streamStage
    for(n in 1:nrow(stageData)){
      stageData$locationID[n] <- gaugeLocData$namedLocation[gaugeLocData$locationStartDate<=stageData$startDate[n]
                                                                &gaugeLocData$locationEndDate>=stageData$startDate[n]
                                                                &gaugeLocData$namedLocation%in%unique(gaugeLocData$namedLocation)]
      stageData$locationID[n] <- gaugeLocData$namedLocation[gaugeLocData$locationStartDate<=stageData$startDate[n]
                                                               &gaugeLocData$locationEndDate>=stageData$startDate[n]
                                                               &gaugeLocData$namedLocation%in%unique(gaugeLocData$namedLocation)]
    }
      
    #Remove any GAG records that are in the DSC data frame and merge
    L1GaugeData <- L1GaugeData[!L1GaugeData$eventID%in%stageData$eventID,]
    L1GaugeData <- merge(L1GaugeData,stageData,all = T)
    
    #Subset out records with no gaugings
    L1GaugeData <- L1GaugeData[!is.na(L1GaugeData$initialStageHeight),]
    
    #Set data types
    L1GaugeData$initialStageHeight <- as.numeric(L1GaugeData$initialStageHeight)
    L1GaugeData$startDate <- as.POSIXct(L1GaugeData$startDate,format=osPubDateFormat)
    L1GaugeData$gaugeHeightOffset <- NA
    L1GaugeData$initialStageHeightWOffset <- NA
    
    # Map L1 gauge data to TROLL data using the HOR
    L1GaugeData$HOR <- NA
    for(h in 1:nrow(L1GaugeData)){
      L1GaugeData$HOR[h] <- gsub("\\..*$","",gaugeLocData$yOffset[gaugeLocData$namedLocation==L1GaugeData$locationID[h]
                                                                  &gaugeLocData$locationStartDate<=format(L1GaugeData$startDate[h],"%Y-%m-%d")
                                                                  &gaugeLocData$locationEndDate+secondsInDay>=format(L1GaugeData$startDate[h],"%Y-%m-%d")])
    }
    
    #Apply offsets for all ranges of gauge data
    namedLocations <- unique(gaugeLocData$namedLocation)
    #Apply offsets for all ranges of gauge data
    for (y in 1:length(namedLocations)) {
      currNamedLocations <- namedLocations[y]
      gaugeLocData_namedLocation <- gaugeLocData[gaugeLocData$namedLocation==currNamedLocations,]
      for(i in 1:length(gaugeLocData_namedLocation$startDate)){
        locStart <- gaugeLocData_namedLocation$locationStartDate[i]
        locEnd <- gaugeLocData_namedLocation$locationEndDate[i]
        dataToApplyOffset <- L1GaugeData$initialStageHeight[L1GaugeData$startDate>=locStart & L1GaugeData$startDate<locEnd]
        if(i==1){
          locOffset <- 0
          L1GaugeData$initialStageHeightWOffset[L1GaugeData$startDate>=locStart & L1GaugeData$startDate<locEnd] <-
            dataToApplyOffset + as.numeric(locOffset)
          L1GaugeData$gaugeHeightOffset[L1GaugeData$startDate>=locStart & L1GaugeData$startDate<locEnd] <-
            as.numeric(locOffset)
        }else{
          locOffset <- gaugeLocData_namedLocation$refElevPlusZ[i] - gaugeLocData_namedLocation$refElevPlusZ[1]
          L1GaugeData$initialStageHeightWOffset[L1GaugeData$startDate>=locStart & L1GaugeData$startDate<locEnd] <-
            dataToApplyOffset + as.numeric(locOffset)
          L1GaugeData$gaugeHeightOffset[L1GaugeData$startDate>=locStart & L1GaugeData$startDate<locEnd] <-
            as.numeric(locOffset)
        }
      }
    }
    
    #Build the empty data frame
    gaugePressureRelationshipData_Names <- c(
      'domainID',
      'siteID',
      'namedLocation',
      'trollLocation',
      'stationHorizontalID',
      'startDate',
      'endDate',
      'calibratedPressMean',
      'calibratedPressObsCount',
      'calibratedPressStdDev',
      'calcWaterColumnHeight',
      'regressionID',
      'sensorStaffGaugeOffset',
      'calculatedStage',
      'gaugeCollectDate',
      'gaugeEventID',
      'gaugeHeight'
    )
    gaugePressureRelationshipData <- data.frame(matrix(data=NA, ncol=length(gaugePressureRelationshipData_Names), nrow=nrow(L1GaugeData)))
    names(gaugePressureRelationshipData) <- gaugePressureRelationshipData_Names
    
    #Order the data and add offset gauge records and collect dates to the L4 table
    L1GaugeData <- L1GaugeData[order(L1GaugeData$startDate),]
    gaugePressureRelationshipData$domainID <- unique(L1GaugeData$domainID)
    gaugePressureRelationshipData$siteID <- site
    gaugePressureRelationshipData$stationHorizontalID <- L1GaugeData$HOR
    gaugePressureRelationshipData$namedLocation <- L1GaugeData$locationID
    gaugePressureRelationshipData$gaugeHeight <-as.numeric(L1GaugeData$initialStageHeightWOffset)
    gaugePressureRelationshipData$gaugeCollectDate <- as.POSIXct(L1GaugeData$startDate,format=osDateFormat)
    gaugePressureRelationshipData$gaugeEventID <- L1GaugeData$eventID
    gaugePressureRelationshipData$startDate <- gaugePressureRelationshipData$gaugeCollectDate
    gaugePressureRelationshipData$endDate <- gaugePressureRelationshipData$gaugeCollectDate
    
    #Pull in L0Pressure for timestamps that match the gauge records
    for (i in 1:nrow(gaugePressureRelationshipData)) {
      #Get a 20 minute mean pressure value for the L4 table
      gaugeCollectDatePlus10Min <- gaugePressureRelationshipData$gaugeCollectDate[i] + secondsIn10Min
      gaugeCollectDateMinus10Min <- gaugePressureRelationshipData$gaugeCollectDate[i] - secondsIn10Min
      pressureDataToAverage <- L0PressureData[L0PressureData$date >= gaugeCollectDateMinus10Min & L0PressureData$date <= (gaugeCollectDatePlus10Min-1),]
      # Skip to next row if there are no records or if all records are within SWE SRFs
      if (all(is.na(pressureDataToAverage$CalPress))|all(pressureDataToAverage$sWatElevFinalQFSciRvw==1&!is.na(pressureDataToAverage$sWatElevFinalQFSciRvw))) {
        print(paste0("No pressure data available at within 20 min of ",gaugePressureRelationshipData$gaugeCollectDate[i],": removing record and moving on to next"))
        next
      }else{
        gaugePressureRelationshipData$calibratedPressMean[i] <- mean(pressureDataToAverage$CalPress, na.rm = T)
        gaugePressureRelationshipData$calibratedPressObsCount[i] <- sum(!is.na(pressureDataToAverage$CalPress))
        gaugePressureRelationshipData$calibratedPressStdDev[i] <- stats::sd(pressureDataToAverage$CalPress, na.rm = T)
      }
    }
    
    #Subset out records with no pressure records
    gaugePressureRelationshipData <- gaugePressureRelationshipData[!is.na(gaugePressureRelationshipData$calibratedPressMean),]
    
    #Convert pressure to water column height above troll
    gaugePressureRelationshipData$waterColumnHeightRaw <- (gaugePressureRelationshipData$calibratedPressMean/(roe * grav)) * convKPatoPa
    
    #Get TROLL location data from the OS system and apply offsets to water column height data
    nlTROLL <- unique(gaugePressureRelationshipData$trollLocation)
    for (y in 1:length(nlTROLL)) {
      currNlTROLL <- nlTROLL[y]
      
      # Get TROLL location data
      trollLocData <- try(geoNEON::getLocBySite(site,type = "AQU",history = T),silent = T)
      # Error handling if the GET was unsuccessful
      if(attr(trollLocData, "class") == "try-error"){
        failureMessage <- "Troll location history could not be retrieved"
        stop(failureMessage)
      }
      trollLocData <- trollLocData[trollLocData$namedLocation%in%unique(L0PressureData$namedLocation),]
      #Error handling if there is no data
      if(nrow(trollLocData)<1){
        failureMessage <- "Zero (0) troll location history records were retrieved"
        stop(failureMessage)
      }
      
      #Set data types
      #Add an endDate to avoid NA and null situations
      if(is.null(trollLocData$locationEndDate)){
        trollLocData$locationEndDate[is.null(trollLocData$locationEndDate)] <- format(as.POSIXct(now),format = osDateFormat)}
      if(any(is.na(trollLocData$locationEndDate))){
        trollLocData$locationEndDate[is.na(trollLocData$locationEndDate)] <- format(as.POSIXct(now),format = osDateFormat)}
      trollLocData$locationStartDate <-  as.POSIXct(trollLocData$locationStartDate,format = osDateFormat)
      trollLocData$locationEndDate <- as.POSIXct(trollLocData$locationEndDate,format = osDateFormat)
      trollLocData$zOffset <- as.numeric(trollLocData$zOffset)
      trollLocData$elevation <- as.numeric(trollLocData$elevation)
      trollLocData <- trollLocData[order(trollLocData$locationEndDate),]
      trollLocData$refElevPlusZ <- trollLocData$elevation + trollLocData$zOffset
      #Apply troll elevation offsets to the calculated water column height in regressionData
      for(n in 1:length(unique(trollLocData$namedLocation))){
        trollLocData_temp <- trollLocData[trollLocData$namedLocation==unique(trollLocData$namedLocation)[n]]
        trollLocData_temp <- trollLocData_temp[order(trollLocData$locationEndDate),]
        for(i in 1:length(trollLocData_temp$locationEndDate)){
          locStart <- trollLocData_temp$locationStartDate[i]
          locEnd <-trollLocData_temp$locationEndDate[i]
          dataToApplyOffset <- gaugePressureRelationshipData$waterColumnHeightRaw[gaugePressureRelationshipData$gaugeCollectDate>=locStart&gaugePressureRelationshipData$gaugeCollectDate<locEnd]
          if(i==1){
            locOffset <- 0
            gaugePressureRelationshipData$calcWaterColumnHeight[gaugePressureRelationshipData$gaugeCollectDate>=locStart&gaugePressureRelationshipData$gaugeCollectDate<locEnd] <- dataToApplyOffset + as.numeric(locOffset)
          }else{
            locOffset <- trollLocData_temp$refElevPlusZ[i] - trollLocData_temp$refElevPlusZ[1]
            gaugePressureRelationshipData$calcWaterColumnHeight[gaugePressureRelationshipData$gaugeCollectDate>=locStart&gaugePressureRelationshipData$gaugeCollectDate<locEnd] <- dataToApplyOffset + as.numeric(locOffset)
          }
        }
      }
    }
    
    #Set data types in the regression data
    regressionData$regressionStartDate <- as.POSIXct(regressionData$regressionStartDate, tz = "UTC", format = osPubDateFormat)
    regressionData$regressionEndDate <- as.POSIXct(regressionData$regressionEndDate, tz = "UTC", format = osPubDateFormat)
    regressionData$regressionSlope <- as.numeric(regressionData$regressionSlope)
    regressionData$regressionIntercept <- as.numeric(regressionData$regressionIntercept)
    
    #Add regression ID to the L0 pressure data
    if(site=="TOOK"){
      regressionData_currLoc <- regressionData[grepl(searchGaugeNamedLocation,regressionData$namedLocation),]
    }else{
      regressionData_currLoc <- regressionData
    }
    for (i in 1:nrow(gaugePressureRelationshipData)) {
      for (j in 1:nrow(regressionData_currLoc)) {
        if (gaugePressureRelationshipData$gaugeCollectDate[i]>=regressionData_currLoc$regressionStartDate[j]&
            gaugePressureRelationshipData$gaugeCollectDate[i]<=regressionData_currLoc$regressionEndDate[j]) {
          gaugePressureRelationshipData$regressionID[i] <- regressionData_currLoc$regressionID[j]
        }
      }
    }
    
    #Calculate stage
    for (i in 1:nrow(gaugePressureRelationshipData)) {
      for (j in 1:nrow(regressionData_currLoc)) {
        if (gaugePressureRelationshipData$regressionID[i]==regressionData_currLoc$regressionID[j]) {
          gaugePressureRelationshipData$calculatedStage[i] <- (gaugePressureRelationshipData$calcWaterColumnHeight[i]*regressionData_currLoc$regressionSlope[j])+regressionData_currLoc$regressionIntercept[j]
        }
      }
    }
    
    #Subset to only those fields in the L4 table
    gaugePressureRelationshipData <- gaugePressureRelationshipData[,gaugePressureRelationshipData_Names]
    
    # For TOOK locations, combine data frames for writing to database.
    if (l==1) {
      gaugePressureRelationshipDataForL4 <- gaugePressureRelationshipData
    }
    if (l>1) {
      gaugePressureRelationshipDataForL4 <- rbind(gaugePressureRelationshipDataForL4,gaugePressureRelationshipData)
    }
  }
  
  #Format dates for transition object
  gaugePressureRelationshipDataForL4$gaugeCollectDate <- format(gaugePressureRelationshipDataForL4$gaugeCollectDate,format=osDateFormat)
  gaugePressureRelationshipDataForL4$startDate <- format(gaugePressureRelationshipDataForL4$startDate,format=osDateFormat)
  gaugePressureRelationshipDataForL4$endDate <- format(gaugePressureRelationshipDataForL4$endDate,format=osDateFormat)
    
  #Format results for transition object to fit into NEON tables
  #Doesn't need curve specific information, just start and end dates that match the water year
  print("Formatting gaugePressureRelationship")
  write.csv(gaugePressureRelationshipDataForL4,paste0(DIRPATH,BAMWS,"gaugePressureRelationship_",site,"_WY",format(WYsearchIntervalEndDate,"%Y"),".csv"),row.names = FALSE)
}
