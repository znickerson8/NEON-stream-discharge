##############################################################################################
#' @title Plot rating curve outputs

#' @author
#' Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr

#' @description This function plots and saves figures from the BaM rating curve prediction
#' outputs

#' @importFrom graphics plot lines polygon points
#' @importFrom grDevices dev.copy2pdf dev.off adjustcolor

#' @param DIRPATH An environment variable that contains the location of the files in
#' the Docker container [string]
#' @param BAMWS An environment variable that contains the location of the BaM config
#' files in the Docker container [string]
#' @param Hgrid A numeric array of the stage values overwhich the rating curve was evaluated
#' and will be plotted [numeric]

#' @return This function writes out configurations and runs BaM in prediction mode

#' @references
#' License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007

#' @export

# changelog and author contributions / copyrights
#   Kaelin M. Cawley (2017-12-07)
#     original creation
##############################################################################################
BaM.RC.out.plot <- function(
  DIRPATH = Sys.getenv("DIRPATH"),
  BAMWS = Sys.getenv("BAMWS"),
  Hgrid
  ){

  library(plotly)

  ### ----- THE CODE BELOW IS FOR A SINGLE RATING CURVE OR THE 1ST SEGMENT OF A RATING CURVE ----- ###

  ## -- Read in output data of the rating curve MCMC predictions
  # Predicted Priors
  Qrc_Prior_spag_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Prior.spag"), header = F)
  Qrc_Prior_env_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Prior.env"), header = T)
  # Predicted Max Post Discharge
  Qrc_Maxpost_spag_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Maxpost.spag"), header = F)
  # Predicted Parametric Uncertainty
  Qrc_ParamU_spag_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_ParamU.spag"), header = F)
  Qrc_ParamU_env_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_ParamU.env"), header = T)
  # Predicted Remnant Uncertainty
  Qrc_TotalU_spag_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_TotalU.spag"), header = F)
  Qrc_TotalU_env_seg1 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_TotalU.env"), header = T)
  # Empirial Gauge and Discharge Pairs
  gaugings_seg1 <- read.table(paste0(DIRPATH, BAMWS, "data/Gaugings.txt"),sep = "\t",header = T)
  gaugings_seg1$Q <- as.numeric(gaugings_seg1$Q) #Convert to cms from lps
  gaugings_seg1$uQ <- as.numeric(gaugings_seg1$uQ) #Convert to cms from lps

  ## -- Get Parameters for plotting
  # Parametric Uncertainty Parameters for Plotting
  pramUForPlottingTop_seg1 <- cbind.data.frame(Hgrid_seg1,Qrc_ParamU_env_seg1$Q_q2.5)
  pramUForPlottingBottom_seg1 <- cbind.data.frame(Hgrid_seg1,Qrc_ParamU_env_seg1$Q_q97.5)
  names(pramUForPlottingTop_seg1) <- c("Hgrid","Q")
  names(pramUForPlottingBottom_seg1) <- c("Hgrid","Q")
  pramUForPlotting_seg1 <- rbind(pramUForPlottingTop_seg1,pramUForPlottingBottom_seg1[dim(pramUForPlottingBottom_seg1)[1]:1,])
  # Remnant Uncertainty Parameters for Plotting
  totalUTop_seg1 <- cbind.data.frame(Hgrid_seg1,Qrc_TotalU_env_seg1$Q_q2.5)
  totalUBottom_seg1 <- cbind.data.frame(Hgrid_seg1,Qrc_TotalU_env_seg1$Q_q97.5)
  names(totalUTop_seg1) <- c("Hgrid","Q")
  names(totalUBottom_seg1) <- c("Hgrid","Q")
  totalUForPlotting_seg1 <- rbind(totalUTop_seg1,totalUBottom_seg1[dim(totalUBottom_seg1)[1]:1,])
  # Prior Uncertainty Parameters for Plotting
  priorTop_seg1 <- cbind.data.frame(Hgrid_seg1, Qrc_Prior_env_seg1$Q_q2.5)
  priorBottom_seg1 <- cbind.data.frame(Hgrid_seg1, Qrc_Prior_env_seg1$Q_q97.5)
  names(priorTop_seg1) <- c("Hgrid","Q")
  names(priorBottom_seg1) <- c("Hgrid","Q")
  priorForPlotting_seg1 <- rbind(priorTop_seg1,priorBottom_seg1[dim(priorBottom_seg1)[1]:1,])

  # ############################################################
  # rcNam <- c('curveID','Hgrid','maxPostQ','pramUTop','pramUBottom','totalUTop','totalUBottom')
  # rcData <- data.frame(matrix(data=NA, ncol=length(rcNam), nrow=181))
  # names(rcData) <- rcNam
  # rcData <- rcData
  # rcData$curveID <- curveID
  # rcData$Hgrid <- Hgrid_seg1
  # rcData$maxPostQ <- Qrc_Maxpost_spag_seg1$V1*1000
  # rcData$pramUTop <- pramUForPlottingTop_seg1$Q*1000
  # rcData$pramUBottom <- pramUForPlottingBottom_seg1$Q*1000
  # rcData$totalUTop <- totalUTop_seg1$Q*1000
  # rcData$totalUBottom <- totalUBottom_seg1$Q*1000
  # # rcData_allWYs <- rcData
  # rcData_allWYs <- rbind(rcData_allWYs,rcData)
  # rcGaugings <- gaugings_seg1
  # rcGaugings$curveID <- curveID
  # # rcGaugings_allWYs <- rcGaugings
  # rcGaugings_allWYs <- rbind(rcGaugings_allWYs,rcGaugings)
  # #
  # rcDataAndGaugings <- list(
  #   rcData_allWYs,
  #   rcGaugings_allWYs
  # )
  # names(rcDataAndGaugings) <- c(
  #   "rcData",
  #   "rcGaugings"
  # )
  # saveRDS(rcDataAndGaugings,"C:/Users/nickerson/Documents/Github/NEON-stream-discharge-divine/L4Discharge/AOSApp/rcPlottingData.rds")
  # #############################################################
  
  ### ----- THE CODE BELOW IS FOR A 2ND SEGMENT OF A RATING CURVE ----- ###

  # ## -- Read in output data of the rating curve MCMC predictions
  # # Predicted Priors
  # Qrc_Prior_spag_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Prior.spag"), header = F)
  # Qrc_Prior_env_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Prior.env"), header = T)
  # # Predicted Max Post Discharge
  # Qrc_Maxpost_spag_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Maxpost.spag"), header = F)
  # # Predicted Parametric Uncertainty
  # Qrc_ParamU_spag_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_ParamU.spag"), header = F)
  # Qrc_ParamU_env_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_ParamU.env"), header = T)
  # # Predicted Remnant Uncertainty
  # Qrc_TotalU_spag_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_TotalU.spag"), header = F)
  # Qrc_TotalU_env_seg2 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_TotalU.env"), header = T)
  # # Empirial Gauge and Discharge Pairs
  # gaugings_seg2 <- read.table(paste0(DIRPATH, BAMWS, "data/Gaugings.txt"),sep = "\t",header = T)
  # gaugings_seg2$Q <- as.numeric(gaugings_seg2$Q) #Convert to cms from lps
  # gaugings_seg2$uQ <- as.numeric(gaugings_seg2$uQ) #Convert to cms from lps
  # 
  # ## -- Get Parameters for plotting
  # # Parametric Uncertainty Parameters for Plotting
  # pramUForPlottingTop_seg2 <- cbind.data.frame(Hgrid_seg2,Qrc_ParamU_env_seg2$Q_q2.5)
  # pramUForPlottingBottom_seg2 <- cbind.data.frame(Hgrid_seg2,Qrc_ParamU_env_seg2$Q_q97.5)
  # names(pramUForPlottingTop_seg2) <- c("Hgrid","Q")
  # names(pramUForPlottingBottom_seg2) <- c("Hgrid","Q")
  # pramUForPlotting_seg2 <- rbind(pramUForPlottingTop_seg2,pramUForPlottingBottom_seg2[dim(pramUForPlottingBottom_seg2)[1]:1,])
  # # Remnant Uncertainty Parameters for Plotting
  # totalUTop_seg2 <- cbind.data.frame(Hgrid_seg2,Qrc_TotalU_env_seg2$Q_q2.5)
  # totalUBottom_seg2 <- cbind.data.frame(Hgrid_seg2,Qrc_TotalU_env_seg2$Q_q97.5)
  # names(totalUTop_seg2) <- c("Hgrid","Q")
  # names(totalUBottom_seg2) <- c("Hgrid","Q")
  # totalUForPlotting_seg2 <- rbind(totalUTop_seg2,totalUBottom_seg2[dim(totalUBottom_seg2)[1]:1,])
  # # Prior Uncertainty Parameters for Plotting
  # priorTop_seg2 <- cbind.data.frame(Hgrid_seg2, Qrc_Prior_env_seg2$Q_q2.5)
  # priorBottom_seg2 <- cbind.data.frame(Hgrid_seg2, Qrc_Prior_env_seg2$Q_q97.5)
  # names(priorTop_seg2) <- c("Hgrid","Q")
  # names(priorBottom_seg2) <- c("Hgrid","Q")
  # priorForPlotting_seg2 <- rbind(priorTop_seg2,priorBottom_seg2[dim(priorBottom_seg2)[1]:1,])

  ### ----- THE CODE BELOW IS FOR A 3RD SEGMENT OF A RATING CURVE ----- ###

  # ## -- Read in output data of the rating curve MCMC predictions
  # # Predicted Priors
  # Qrc_Prior_spag_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Prior.spag"), header = F)
  # Qrc_Prior_env_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Prior.env"), header = T)
  # # Predicted Max Post Discharge
  # Qrc_Maxpost_spag_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_Maxpost.spag"), header = F)
  # # Predicted Parametric Uncertainty
  # Qrc_ParamU_spag_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_ParamU.spag"), header = F)
  # Qrc_ParamU_env_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_ParamU.env"), header = T)
  # # Predicted Remnant Uncertainty
  # Qrc_TotalU_spag_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_TotalU.spag"), header = F)
  # Qrc_TotalU_env_seg3 <- read.table(paste0(DIRPATH, BAMWS, "Qrc_TotalU.env"), header = T)
  # # Empirial Gauge and Discharge Pairs
  # gaugings_seg3 <- read.table(paste0(DIRPATH, BAMWS, "data/Gaugings.txt"),sep = "\t",header = T)
  # gaugings_seg3$Q <- as.numeric(gaugings_seg3$Q) #Convert to cms from lps
  # gaugings_seg3$uQ <- as.numeric(gaugings_seg3$uQ) #Convert to cms from lps
  #
  # ## -- Get Parameters for plotting
  # # Parametric Uncertainty Parameters for Plotting
  # pramUForPlottingTop_seg3 <- cbind.data.frame(Hgrid_seg3,Qrc_ParamU_env_seg3$Q_q2.5)
  # pramUForPlottingBottom_seg3 <- cbind.data.frame(Hgrid_seg3,Qrc_ParamU_env_seg3$Q_q97.5)
  # names(pramUForPlottingTop_seg3) <- c("Hgrid","Q")
  # names(pramUForPlottingBottom_seg3) <- c("Hgrid","Q")
  # pramUForPlotting_seg3 <- rbind(pramUForPlottingTop_seg3,pramUForPlottingBottom_seg3[dim(pramUForPlottingBottom_seg3)[1]:1,])
  # # Remnant Uncertainty Parameters for Plotting
  # totalUTop_seg3 <- cbind.data.frame(Hgrid_seg3,Qrc_TotalU_env_seg3$Q_q2.5)
  # totalUBottom_seg3 <- cbind.data.frame(Hgrid_seg3,Qrc_TotalU_env_seg3$Q_q97.5)
  # names(totalUTop_seg3) <- c("Hgrid","Q")
  # names(totalUBottom_seg3) <- c("Hgrid","Q")
  # totalUForPlotting_seg3 <- rbind(totalUTop_seg3,totalUBottom_seg3[dim(totalUBottom_seg3)[1]:1,])
  # # Prior Uncertainty Parameters for Plotting
  # priorTop_seg3 <- cbind.data.frame(Hgrid_seg3, Qrc_Prior_env_seg3$Q_q2.5)
  # priorBottom_seg3 <- cbind.data.frame(Hgrid_seg3, Qrc_Prior_env_seg3$Q_q97.5)
  # names(priorTop_seg3) <- c("Hgrid","Q")
  # names(priorBottom_seg3) <- c("Hgrid","Q")
  # priorForPlotting_seg3 <- rbind(priorTop_seg3,priorBottom_seg3[dim(priorBottom_seg3)[1]:1,])


  ### ----- Plot the rating curve with uncertainty on a linear scale ----- ###

  # ## -- Base R -- ##
  # # Set default plotting margins
  # test.par <- par(mfrow=c(1,1))
  # par(mar=c(5.1,4.1,4.1,2.1))
  # plot(Hgrid,
  #      Qrc_Maxpost_spag$V1,
  #      type = "l",
  #      xlim = c(min(Hgrid),max(Hgrid)),
  #      #xlim = WY2019xlim,
  #      #ylim = c(min(priorForPlotting$Q),max(priorForPlotting$Q)), #Cannot plot due to impossible spaghettis
  #      ylim = c(min(totalUForPlotting$Q),max(totalUForPlotting$Q)), #Total uncertainty replacing priors as the y-axis limits
  #      #ylim = WY2019ylim,
  #      xlab = "Stage (m)",
  #      ylab = "Discharge (cms)")
  # # Work from background to foreground to add layers to the plot
  # # Red shades for uncertainty of a rating #1 for WY
  # polygon(totalUForPlotting$Hgrid,totalUForPlotting$Q, col = "red", border = NA)
  # polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = adjustcolor("lightpink",alpha.f = 0.7), border = NA)
  # # # Blue shades for uncertainty of a rating #2 for WY
  # # polygon(totalUForPlotting$Hgrid,totalUForPlotting$Q, col = "blue", border = NA)
  # # polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = adjustcolor("lightblue",alpha.f = 0.7), border = NA)
  # # # Green shades for uncertainty of a rating #3 for WY
  # # polygon(totalUForPlotting$Hgrid,totalUForPlotting$Q, col = "green4", border = NA)
  # # polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = adjustcolor("green",alpha.f = 0.7), border = NA)
  # # Add lines and points for rating curve and gauge/discharge pairs, respectively
  # lines(Hgrid,Qrc_Maxpost_spag$V1, col = "black", lwd = 2)
  # points(gaugings$H,gaugings$Q, pch = 19, col = "black")
  # # Place min/max calculated stage in the plot
  # abline(v=minCalcStage,lwd=3,lty=2)
  # abline(v=maxCalcStage,lwd=3,lty=2)
  # # Save out the plot to the destination file path
  # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"priorAndPostRatingCurves_linearScale_MAYF2018_updatedControls3_minMacCalcStageBounds.pdf",sep = "/"), width = 16, height = 9)
  # dev.off()

  ## -- plotly -- ##
  setwd(paste0(DIRPATH,BAMWS))
  (plotly <- plot_ly()%>%

      ### ----- USE THE BELOW AESTHETICS FOR A SINGLE RATING CURVE OR THE 1ST SEGMENT OF A RATING CURVE ------ ###
      # Total Uncertainty
      add_trace(x=totalUTop_seg1$Hgrid,y=totalUTop_seg1$Q*1000,name='Remn U Top 1',type='scatter',mode='line',line=list(color='red'),legendgroup='group1')%>%
      add_trace(x=totalUBottom_seg1$Hgrid,y=totalUBottom_seg1$Q*1000,name='Remn U Bottom 1',type='scatter',mode='line',fill='tonexty',fillcolor='red',line=list(color='red'),legendgroup='group1')%>%
      # Parametric Uncertainty
      add_trace(x=pramUForPlottingTop_seg1$Hgrid,y=pramUForPlottingTop_seg1$Q*1000,name='Para U Top 1',type='scatter',mode='line',line=list(color='lightpink'),legendgroup='group1')%>%
      add_trace(x=pramUForPlottingBottom_seg1$Hgrid,y=pramUForPlottingBottom_seg1$Q*1000,name='Para U Bottom 1',type='scatter',mode='line',fill='tonexty',fillcolor='lightpink',line=list(color='lightpink'),legendgroup='group1')%>%
      # Max Post Q
      add_trace(x=Hgrid_seg1,y=Qrc_Maxpost_spag_seg1$V1*1000,name='Max Post Q 1',type='scatter',mode='line',line=list(color='black'),legendgroup='group2')%>%
      # Empirical H/Q Pairs
      add_trace(x=gaugings_seg1$H,y=gaugings_seg1$Q,name='Empirial H/Q Pairs 1',type='scatter',mode='markers',marker=list(color='black'),legendgroup='group2')%>%

      # ## ----- USE THE BELOW AESTHETICS FOR A 2ND SEGMENT OF A RATING CURVE ------ ###
      # # Total Uncertainty
      # add_trace(x=totalUTop_seg2$Hgrid,y=totalUTop_seg2$Q*1000,name='Remn U Top 2',type='scatter',mode='line',line=list(color='blue'),legendgroup='group3')%>%
      # add_trace(x=totalUBottom_seg2$Hgrid,y=totalUBottom_seg2$Q*1000,name='Remn U Bottom 2',type='scatter',mode='line',fill='tonexty',fillcolor='blue',line=list(color='blue'),legendgroup='group3')%>%
      # # Parametric Uncertainty
      # add_trace(x=pramUForPlottingTop_seg2$Hgrid,y=pramUForPlottingTop_seg2$Q*1000,name='Para U Top 2',type='scatter',mode='line',line=list(color='lightblue'),legendgroup='group3')%>%
      # add_trace(x=pramUForPlottingBottom_seg2$Hgrid,y=pramUForPlottingBottom_seg2$Q*1000,name='Para U Bottom 2',type='scatter',mode='line',fill='tonexty',fillcolor='lightblue',line=list(color='lightblue'),legendgroup='group3')%>%
      # # Max Post Q
      # add_trace(x=Hgrid_seg2,y=Qrc_Maxpost_spag_seg2$V1*1000,name='Max Post Q 2',type='scatter',mode='line',line=list(color='black'),legendgroup='group4')%>%
      # # Empirical H/Q Pairs
      # add_trace(x=gaugings_seg2$H,y=gaugings_seg2$Q,name='Empirial H/Q Pairs 2',type='scatter',mode='markers',marker=list(color='black'),legendgroup='group4')%>%

      # ## ----- USE THE BELOW AESTHETICS FOR A 3RD SEGMENT OF A RATING CURVE ------ ###
      # # Total Uncertainty
      # add_trace(x=totalUTop_seg3$Hgrid,y=totalUTop_seg3$Q,name='Remn U Top 3',type='scatter',mode='line',line=list(color='green'),legendgroup='group5')%>%
      # add_trace(x=totalUBottom_seg3$Hgrid,y=totalUBottom_seg3$Q,name='Remn U Bottom 3',type='scatter',mode='line',fill='tonexty',fillcolor='green',line=list(color='green'),legendgroup='group5')%>%
      # # Parametric Uncertainty
      # add_trace(x=pramUForPlottingTop_seg3$Hgrid,y=pramUForPlottingTop_seg3$Q,name='Para U Top 3',type='scatter',mode='line',line=list(color='lightgreen'),legendgroup='group5')%>%
      # add_trace(x=pramUForPlottingBottom_seg3$Hgrid,y=pramUForPlottingBottom_seg3$Q,name='Para U Bottom 3',type='scatter',mode='line',fill='tonexty',fillcolor='lightgreen',line=list(color='lightgreen'),legendgroup='group5')%>%
      # # Max Post Q
      # add_trace(x=Hgrid_seg3,y=Qrc_Maxpost_spag_seg3$V1,name='Max Post Q 3',type='scatter',mode='line',line=list(color='black'),legendgroup='group6')%>%
      # # Empirical H/Q Pairs
      # add_trace(x=gaugings_seg3$H,y=gaugings_seg3$Q,name='Empirial H/Q Pairs 3',type='scatter',mode='markers',marker=list(color='black'),legendgroup='group6')%>%

      ### ----- BELOW ARE AESTHETICS FOR ALL PLOTS ----- ###
      # # Min/Max Calc H
      # add_segments(x=minCalcH,xend=minCalcH,y=0,yend=max(totalUForPlotting_seg1$Q)*1000,name='Min Calc Q',showlegend=F,line=list(color='black',dash='dash'))%>%
      # add_segments(x=maxCalcH,xend=maxCalcH,y=0,yend=max(totalUForPlotting_seg1$Q)*1000,name='Max Calc Q',showlegend=F,line=list(color='black',dash='dash'))%>%
      #Plot Layout
      layout(title=paste0(gsub("-[1-3]{1}","",unique(sampledParameters$curveID))," -- Rating Curve w/ Uncertainty -- Linear"),
        xaxis=list(title="Stage (m)"),
        yaxis=list(title="Discharge (cms)"),
        updatemenus=list(
          list(
            type='buttons',
            buttons=list(
              list(label='linear',
                method='relayout',
                args=list(list(yaxis=list(type='linear')))),
              list(label='log',
                method='relayout',
                args=list(list(yaxis=list(type='log')))))))))
  #Save out the plotly image
  htmlwidgets::saveWidget(as_widget(plotly), paste0(gsub("-[1-3]{1}","",unique(sampledParameters$curveID)),"_posteriorRC.html"))

  # # #Plot rating curve with gaugings
  # # plot(Hgrid,
  # #      Qrc_Maxpost_spag$V1,
  # #      type = "l",
  # #      xlim = c(min(Hgrid),max(Hgrid)),
  # #      ylim = c(min(pramUForPlotting$Q,gaugings$H),max(pramUForPlotting$Q,gaugings$H)),
  # #      xlab = "Stage (m)",
  # #      ylab = "Discharge (cms)")
  # # polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = "lightpink", border = NA)
  # # lines(Hgrid,Qrc_Maxpost_spag$V1, col = "red", lwd = 2)
  # # points(gaugings$H,gaugings$Q, col="red", pch = 19)
  # #
  # # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"ratingCurveWithGaugings_linearScale_WALK2018.pdf",sep = "/"), width = 16, height = 9)
  # # dev.off()
  #
  # #Plot log-scale
  # test.par <- par(mfrow=c(1,1))
  # par(mar=c(5.1,4.1,4.1,2.1))
  # Qrc_Maxpost_spag$V1[Qrc_Maxpost_spag$V1 <= 0] <- 0.000000001
  # plot(Hgrid,
  #      Qrc_Maxpost_spag$V1,
  #      type = "l",
  #      log = "xy",
  #      xlim = c(9.708,max(Hgrid)),
  #      #xlim = c(0.008,max(Hgrid)),
  #      ylim = c(0.0005,max(totalUForPlotting$Q)), #Total uncertainty replacing priors as the y-axis limits
  #      xlab = "Stage (m)",
  #      ylab = "Discharge (cms)")
  #
  # #Work from background to foreground
  # totalUForPlotting$Qlog <- totalUForPlotting$Q
  # totalUForPlotting$Qlog[totalUForPlotting$Q<=0] <- 0.000001
  # #polygon(priorForPlotting$Hgrid, priorForPlotting$Q, col = "royalblue1", border = NA) #Cannot plot priors due to impossible spaghettis
  # polygon(totalUForPlotting$Hgrid,totalUForPlotting$Qlog, col = "red", border = NA)
  # polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = adjustcolor("lightpink",alpha.f = 0.7), border = NA)
  # #polygon(totalUForPlotting$Hgrid,totalUForPlotting$Qlog, col = "blue", border = NA)
  # #polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = adjustcolor("lightblue",alpha.f = 0.7), border = NA)
  # #polygon(totalUForPlotting$Hgrid,totalUForPlotting$Qlog, col = "green4", border = NA)
  # #polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = adjustcolor("green",alpha.f = 0.7), border = NA)
  #
  # #lines(Hgrid,Qrc_Prior_env$Q_Median, col = "blue", lwd = 2) #Cannot plot priors due to impossible spaghettis
  # lines(Hgrid,Qrc_Maxpost_spag$V1, col = "black", lwd = 2)
  # points(gaugings$H+10,gaugings$Q, pch = 19, col = "black")
  # #lines(Hgrid,Qrc_Maxpost_spag$V1, col = "blue", lwd = 2)
  # #points(gaugings$H,gaugings$Q, pch = 19, col = "blue")
  # #lines(Hgrid,Qrc_Maxpost_spag$V1, col = "red", lwd = 2)
  # #points(gaugings$H,gaugings$Q, pch = 19, col = "red")
  #
  # # Place min/max calculated stage in the plot
  # abline(v=(minCalcStage+10),lwd=3,lty=2)
  # abline(v=(maxCalcStage+10),lwd=3,lty=2)
  #
  # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"priorAndPostRatingCurves_logScale_MAYF2018_updatedControls3_minMaxCalcStageBounds.pdf",sep = "/"), width = 16, height = 9)
  # dev.off()
  #
  # #Plot rating curve with gaugings in log scale
  # plot(Hgrid,
  #      Qrc_Maxpost_spag$V1,
  #      type = "l",
  #      log = "xy",
  #      xlim = c(0.008,max(Hgrid)),
  #      #ylim = c(0.0005,max(priorForPlotting$Q)), #Cannot plot due to impossible spaghettis
  #      ylim = c(0.0005,max(totalUForPlotting$Q)), #Total uncertainty replacing priors as the y-axis limits
  #      xlab = "Stage (m)",
  #      ylab = "Discharge (cms)")
  # polygon(pramUForPlotting$Hgrid,pramUForPlotting$Q, col = "lightpink", border = NA)
  # lines(Hgrid,Qrc_Maxpost_spag$V1, col = "red", lwd = 2)
  # points(gaugings$H,gaugings$Q, col = "red",pch = 19)
  # #lines(Hgrid,Qrc_Maxpost_spag$V1, col = "blue", lwd = 2)
  # #points(gaugings$H,gaugings$Q, col = "blue",pch = 19)
  # #lines(Hgrid,Qrc_Maxpost_spag$V1, col = "green4", lwd = 2)
  # #points(gaugings$H,gaugings$Q, col = "green4",pch = 19)
  #
  # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"ratingCurveWithGaugings_logScale_MAYF2019_updatedControls_minMaxCalcH.pdf",sep = "/"), width = 16, height = 9)
  # dev.off()
  #
  # #Plot spaghettis
  # #Total Uncertainty quite large
  # plot(Hgrid,
  #      Qrc_TotalU_spag[,1],
  #      type = "l",
  #      col = "red",
  #      ylim = c(min(Qrc_TotalU_spag),max(Qrc_TotalU_spag)),
  #      ylab = "Spaghettis with total uncertainty")
  # for(i in 2:ncol(Qrc_TotalU_spag)){
  #   lines(Hgrid,Qrc_TotalU_spag[,i], col = "red")
  # }
  # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"spaghettisWithTotalU.pdf",sep = "/"), width = 16, height = 9)
  # dev.off()
  # #Parametric Uncertainty
  # plot(Hgrid,
  #      Qrc_ParamU_spag[,1],
  #      type = "l",
  #      col = "lightpink",
  #      ylim = c(min(Qrc_ParamU_spag),max(Qrc_ParamU_spag)),
  #      ylab = "Spaghettis with parametric uncertainty")
  # for(i in 2:ncol(Qrc_ParamU_spag)){
  #   lines(Hgrid,Qrc_ParamU_spag[,i], col = "lightpink")
  # }
  # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"spaghettisWithParamU.pdf",sep = "/"), width = 16, height = 9)
  # dev.off()
  # #Prior spaghettis
  # plot(Hgrid,
  #      Qrc_Prior_spag[,1],
  #      type = "l",
  #      col = "royalblue1",
  #      ylim = c(min(Qrc_Prior_spag),max(Qrc_Prior_spag)), #Cannot plot due to impossible spaghettis
  #      #ylim = c(-1,35), #Plotted with arbitrary limits
  #      ylab = "Prior spaghettis")
  # for(i in 2:ncol(Qrc_Prior_spag)){
  #   lines(Hgrid,Qrc_Prior_spag[,i], col = "royalblue1")
  # }
  # dev.copy2pdf(file = paste(DIRPATH,BAMWS,"priorSpaghettis.pdf",sep = "/"), width = 16, height = 9)
  # dev.off()

}
