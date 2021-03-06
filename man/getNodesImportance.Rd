% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drivingforce.R
\name{getNodesImportance}
\alias{getNodesImportance}
\title{Calculate the importance of a node in a graph}
\usage{
getNodesImportance(graph, nodes = NULL, measurements = c("DFC", "DCC",
  "DDC", "DC", "CC", "BC"), componentType = "weak", pathType = "all",
  directed = FALSE, normalized = TRUE)
}
\arguments{
\item{componentType}{character string, either weak or strong; for directed graphs, weak implies weakly, strong strongly connected components to search.}

\item{pathType}{defined the types of the paths used for measuring the distance in directed graphs.}

\item{directed}{whether the graph is directed or undirected;}

\item{normalized}{whether to rescale the results to [0,1]}

\item{g}{an igraph object}

\item{measurement}{DFC, DCC, DDC, DC, CC, BC}
}
\value{
a importance matrix
}
\description{
Calculate the importance of a node in a graph
}
