% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/overlaps.R
\name{count_gr_overlaps}
\alias{count_gr_overlaps}
\title{Count overlaps between two GRanges objects.}
\usage{
count_gr_overlaps(fragments, target_GRanges, binarize, aggregate, sparse)
}
\arguments{
\item{target_GRanges}{A single GenomicRanges object to use as a set of target regions}

\item{binarize}{A logical object indicating whether or not to binarize overlaps (count any number of overlaps between query and a given target region as 1). Default is TRUE.}

\item{aggregate}{A logical object indicating whether the results should be aggregated to a vector with the sum of counts for each query to all target regions. Default is FALSE.}

\item{sparse}{A logical object indicating whether the results should be a sparse matrix (TRUE) or a full matrix (FALSE). Default is TRUE.}

\item{query_fragments}{A GenomicRanges object or a list of GenomicRanges objects to use as a query}
}
\value{
If aggregate = TRUE, a vector of counts. If aggregate is FALSE and sparse is FALSE, a vector of values. If sparse is true, a list to use to build a dgCMatrix object.
}
\description{
Count overlaps between two GRanges objects.
}
