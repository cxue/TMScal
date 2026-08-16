\name{TMScal-package}
\alias{TMScal-package}
\alias{TMScal}
\docType{package}
\title{Tumor Mutational Spectrum Prognostic Calculator}
\description{
  A tool for predicting prognostic risks using Tumor Mutational Spectrum (TMS) 
  from whole exome sequencing (WES) data of MSS (microsatellite stable) 
  cancer patients.
}
\details{
  The method uses 192 trinucleotide mutational pathways with adaptive 
  stability-based feature selection. Provides functions for:
  \itemize{
    \item Data preparation and format conversion
    \item Pathway stability analysis
    \item Feature selection (stability + Cox regression)
    \item Model training and validation
    \item Prognostic risk prediction
    \item Visualization (KM curves, stability plots)
  }
}
\author{Xue Cheng}
\keyword{package}