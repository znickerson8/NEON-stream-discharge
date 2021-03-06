% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/frmt.gaug.disc.mea.file.R
\name{frmt.gaug.disc.mea.file}
\alias{frmt.gaug.disc.mea.file}
\title{Format Stage-Discharge File}
\usage{
frmt.gaug.disc.mea.file(dataFrame, metadata, curveIDData)
}
\arguments{
\item{dataFrame}{A dataframe containing data to be formatted [dataframe]}

\item{metadata}{A list containing the transition metadata that includes: domain,site,
startDateFormatted,endDateFormatted,namedLocationName,numCtrls,numCurves,waterYear [list]}

\item{curveIDData}{A dataframe containing rating curve-specific start and end datetimes [dataframe]}
}
\value{
This function returns a dataframe formatted identical to the similar publication
table in the Stage-discharge rating curve (DP4.00133.001) data product.
}
\description{
This function takes a dataframe and .
}
\references{
License: GNU AFFERO GENERAL PUBLIC LICENSE Version 3, 19 November 2007
}
\author{
Kaelin M. Cawley \email{kcawley@battelleecology.org} \cr
Zachary L. Nickerson \email{nickerson@battelleecology.org} \cr
}
