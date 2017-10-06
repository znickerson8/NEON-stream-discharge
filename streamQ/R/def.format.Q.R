##############################################################################################
#' @title Stream discharge formatting

#' @author 
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function reads in and formats data from the NEON Salt-based Discharge 
#' data product to calculate stream discharge from slug and constant-rate injection. 
#' No need to unzip the downloaded files, just place them all in the smae directory.

#' @importFrom neonDataStackR stackByTable
#' @importFrom utils read.csv

#' @param dataDir User identifies the directory that contains the zipped data

#' @return This function returns one data frame formatted for use with def.calc.Q.slug.R 
#' and def.calc.Q.constantRate.R to calculate stream discharge

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @keywords surface water, streams, rivers, discharge, tracer, salt-based, slug, constant rate

#' @examples
#' #where the data .zip file is in the working directory and has the default name, 
#' sbdFormatted <- def.format.Q()
#' #where the data.zip file is in the downloads folder and has default name, 
#' sbdFormatted <- def.format.Q(dataDir = path.expand("~/NEON_discharge-stream-saltbased.zip"))
#' #where the data.zip file is in the downloads folder and has a specified name,
#' sbdFormatted <- def.format.Q(dataDir = path.expand("~/Downloads/non-standard-name.zip"))
#' #Using the example data in this package
#' dataDirectory <- paste(path.package("neonStreamQ"),"inst\\extdata", sep = "\\")
#' sbdFormatted <- def.format.Q(dataDir = dataDirectory)

#' @seealso def.calc.Q.slug.R and def.calc.Q.inj.R to calculate stream discharge

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-08-03)
#     original creation
##############################################################################################
#This code is for calculating salt-based discharge for a slug
def.format.Q <- function(
  dataDir = paste0(getwd(),"/NEON_discharge-stream-saltbased.zip")
) {
  
  
  #Stack field and external lab data
  if(!dir.exists(substr(dataDir, 1, (nchar(dataDir)-4)))){
    stackByTable(dataDir)
  }
  
  #Read in stacked data
  backgroundDataLogger <- read.csv(
    paste(gsub("\\.zip","",dataDir), "stackedFiles", "sbd_backgroundFieldCondData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  backgroundDataSalt <- read.csv(
    paste(gsub("\\.zip","",dataDir), "stackedFiles", "sbd_backgroundFieldSaltData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  fieldDataSite <- read.csv(
    paste(gsub("\\.zip","",dataDir), "stackedFiles", "sbd_fieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  fieldDataSite$namedLocation <- NULL #So that the merge goes smoothly
  
  plateauDataCond <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","sbd_plateauMeasurementFieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  plateauDataSalt <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","sbd_plateauSampleFieldData.csv", sep = "/"), 
    stringsAsFactors = F)
  
  externalLabDataSalt <- read.csv(
    paste(gsub("\\.zip","",dataDir),"stackedFiles","sbd_externalLabDataSalt.csv", sep = "/"), 
    stringsAsFactors = F)
  
  #Remove two bad duplicates
  externalLabDataSalt <- externalLabDataSalt[!(externalLabDataSalt$saltSampleID == 'POSE.00.20140811.TCR' & externalLabDataSalt$laboratoryName == 'Loeke Lab at University of Kansas'),]
  externalLabDataSalt <- externalLabDataSalt[!(externalLabDataSalt$saltSampleID == 'POSE.00.20141007.TCR' & externalLabDataSalt$laboratoryName == 'Loeke Lab at University of Kansas'),]
  
  #A little error handling until there is a de-duping function
  extDuplicates <- externalLabDataSalt$saltSampleID[duplicated(externalLabDataSalt$saltSampleID)]
  if(length(extDuplicates) > 0){
    stop(paste("Reconcile duplicate external lab sample(s):",paste0(extDuplicates, collapse = ", ")))
  }
  
  #Merge the backgroundDataLogger and fieldDataSite tables
  loggerSiteData <- merge(backgroundDataLogger, 
                          fieldDataSite, 
                          by = c('siteID', 'startDate'), 
                          all = T)
  
  #Create input file for Q calculations
  outputDFNames <- c(
    'siteID',
    'namedLocation',
    'startDate',
    'injectionType',
    
    'slugTracerMass',
    'slugPourTime',
    
    'injectateSampleID',
    'injConc',
    'constantRateTracerMass',
    'carboyVolume',
    'dripStartTime',
    'injRate',
    
    'backgroundConc',
    'plateauConc',
    'hoboSampleID'
  )
  outputDF <- data.frame(matrix(data=NA, ncol=length(outputDFNames), nrow=length(loggerSiteData$siteID)))
  names(outputDF) <- outputDFNames
  
  #Set default QFs
  outputDF$injConcCalcQF <- 0
  
  #Fill in the fields from the loggerSiteData table
  for(i in seq(along = names(outputDF))){
    if(names(outputDF)[i] %in% names(loggerSiteData)){
      outputDF[,i] <- loggerSiteData[,which(names(loggerSiteData) == names(outputDF)[i])]
    }
  }
  
  #Use dripStartTime where slugPoutTime is missing
  outputDF$slugPourTime[is.na(outputDF$slugPourTime) & outputDF$injectionType != 'NaCl' & !is.na(outputDF$dripStartTime)] <-
    outputDF$dripStartTime[is.na(outputDF$slugPourTime) & outputDF$injectionType != 'NaCl' & !is.na(outputDF$dripStartTime)]
  
  #Calculate mean injection rate
  outputDF$injRate <- apply(loggerSiteData[,which(names(loggerSiteData) == "dripRateStart"|names(loggerSiteData) == "dripRateEnd")], 1, mean, na.rm = T)
  
  #Fill in fields that are means of reps
  plateauFields <- c('specificConductanceRep1',
                     'specificConductanceRep2',
                     'specificConductanceRep3',
                     'specificConductanceRep4',
                     'specificConductanceRep5')
  
  for(i in seq(along = outputDF$siteID)){
    
    siteID <- outputDF$siteID[i]
    startDate <- outputDF$startDate[i]
    station <- outputDF$namedLocation[i]
    stationType <- substr(station, 6, nchar(station))
    
    repRegex <- switch(stationType,
           "AOS.reaeration.station.01" = "\\.0[12345]\\.",
           "AOS.reaeration.station.02" = "\\.0[6789]\\.|\\.10\\.",
           "AOS.reaeration.station.03" = "\\.1[12345]\\.",
           "AOS.reaeration.station.04" = "\\.1[6789]\\.|\\.20\\."
    )
    
    #Use last background startDate where slugPourTime is missing
    if(is.na(outputDF$slugPourTime[i])){
      outputDF$slugPourTime[i] <- min(backgroundDataSalt$startDate[backgroundDataSalt$siteID == siteID &
                                                                     backgroundDataSalt$startDate == startDate])
    }
    
    #Populate injectate concentration
    try(outputDF$injConc[i] <- unique(externalLabDataSalt$finalConcentration[
      externalLabDataSalt$saltSampleID == outputDF$injectateSampleID[i]]), silent = T)
    
    #Find plateau values to average
    pValues <- c(t(plateauDataCond[plateauDataCond$startDate == startDate &
                                     plateauDataCond$namedLocation == station,
                                   names(plateauDataCond) %in% plateauFields]))
    #Remove outliers TBD
    #Calculate the mean
    outputDF$plateauCond[i] <- mean(pValues, na.rm = T)
    
    #Fill in background concentration data
    try(outputDF$backgroundConc[i] <- 
        unique(externalLabDataSalt$finalConcentration[
            externalLabDataSalt$namedLocation == station &
              externalLabDataSalt$startDate == startDate & 
              grepl(paste0(".B",substr(station,nchar(station),nchar(station)),"."),
                    externalLabDataSalt$saltSampleID)]), silent = T)
    
    #Fill in plateau concentration data for constant rate injection
    pSaltConc <- externalLabDataSalt$finalConcentration[
      externalLabDataSalt$namedLocation == station &
        externalLabDataSalt$startDate == startDate & 
        grepl(repRegex, externalLabDataSalt$saltSampleID)]
    
    #Remove outliers TBD
    #Calculate the mean plateau concentration
    outputDF$plateauConc[i] <- mean(pSaltConc, na.rm = T)
  }
  
  #Flag calculated values
  outputDF$injConcCalcQF[is.na(outputDF$injConc) & 
                           outputDF$injectionType != 'model' & 
                           !is.na(outputDF$constantRateTracerMass) & 
                           !is.na(outputDF$carboyVolume)] <- 1
  
  #Calculate concentration of constant-rate injectate if no external lab data
  outputDF$injConc[is.na(outputDF$injConc) & 
                     outputDF$injectionType != 'model' & 
                     !is.na(outputDF$constantRateTracerMass) & 
                     !is.na(outputDF$carboyVolume)] <-
    outputDF$constantRateTracerMass[is.na(outputDF$injConc) & 
                                      outputDF$injectionType != 'model' & 
                                      !is.na(outputDF$constantRateTracerMass) & 
                                      !is.na(outputDF$carboyVolume)] * 1000 / 
    outputDF$carboyVolume[is.na(outputDF$injConc) & 
                            outputDF$injectionType != 'model' & 
                            !is.na(outputDF$constantRateTracerMass) & 
                            !is.na(outputDF$carboyVolume)]
  

  return(outputDF)
}
                                           