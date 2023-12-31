% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clustering.R
\name{atac_clustering}
\alias{atac_clustering}
\title{Perform Jaccard-Louvain clustering of ATAC LSI results}
\usage{
atac_clustering(lsi_mat, k = 30, mutual = FALSE, radius = NULL)
}
\arguments{
\item{lsi_mat}{A matrix of LSI results, as generated by atac_lsi()}

\item{k}{A numeric object specifying the number of K-nearest neighbors to use for clustering}

\item{mutual}{A logical object specifying whether or not to exlusively use mutual nearest neighbors for clustering}

\item{radius}{A maximum radius to use to restrict neighbor distance. Default is NULL, which will not impose this restriction.}
}
\value{
a list object with a cluster factor object (cl), modularity scores for each cluster, and the results generated by igraph::cluster_louvain()
}
\description{
Perform Jaccard-Louvain clustering of ATAC LSI results
}
