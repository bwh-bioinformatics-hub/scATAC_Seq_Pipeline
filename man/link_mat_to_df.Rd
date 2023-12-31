% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/archr.R
\name{link_mat_to_df}
\alias{link_mat_to_df}
\title{Convert a sparse linkage matrix to an edge data.frame for graph analysis}
\usage{
link_mat_to_df(link_mat)
}
\arguments{
\item{link_mat}{a matrix with sources as columns, targets as rows, and weights as values.}
}
\value{
a data.frame with columns "from", "to", and "weight".
}
\description{
Convert a sparse linkage matrix to an edge data.frame for graph analysis
}
