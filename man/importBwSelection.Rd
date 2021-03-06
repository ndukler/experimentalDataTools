% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bigWigLib.R
\name{importBwSelection}
\alias{importBwSelection}
\title{Imports bigWig(s)}
\usage{
importBwSelection(expDes, gReg.gr, as.type = "RleList", nthreads = 1)
}
\arguments{
\item{expDes}{experimentalDesign object}

\item{gReg.gr}{either a GRangesList or GRanges object containing genomic coordinates to import}

\item{as.type}{how to return the imported data}

\item{nthreads}{number of cores to use to import data (unimplemented)}
}
\description{
Import bigwig files for later data retrieval.
}
