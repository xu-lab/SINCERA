% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drivingforce.R
\name{drivingforce.getLCC}
\alias{drivingforce.getLCC}
\title{Extract largest connected component (LCC) of the cell group specific TRN}
\usage{
drivingforce.getLCC(object, groups = NULL, thresh = 0.05, ...)
}
\arguments{
\item{object}{A sincera object}

\item{groups}{The cell groups for LCC extraction}

\item{thresh}{The threshold for significant interactions}
}
\value{
An updated sincera object, use getDF with name="lcc" to access the lcc
}
\description{
Extract largest connected component (LCC) of the cell group specific TRN
}
