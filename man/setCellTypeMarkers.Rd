% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sincera.R
\name{setCellTypeMarkers}
\alias{setCellTypeMarkers}
\title{Add cell type marker information into sincera}
\usage{
setCellTypeMarkers(object, value)
}
\arguments{
\item{object}{A sincera object}

\item{value}{A data frame containing the markers and their corresponding cell types. The first column is cell type, and the second column is marker gene name.}
}
\description{
Add cell type marker information into sincera
}
\details{
The value is a data frame with the first column containing the cell types and the second column containing the markers, e.g.,
 TYPE   SYMBOL
  AT2   SFTPB
  AT2   ABCA3
  AT2   SLC34A2
  AT2   LPCAT1
 Basal  KRT5
 Basal  KRT14
}
