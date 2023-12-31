% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{convert_fragments_bed}
\alias{convert_fragments_bed}
\title{Convert a list of fragment GenomicRanges to data.table}
\usage{
convert_fragments_bed(fragments, n_threads = 1)
}
\arguments{
\item{fragments}{A list of GenomicRanges generated by convert_fragments_gr()}

\item{n_threads}{A numeric value indicating how many threads to use for parallelization. Only works on UNIX-like systems. Default is 1.}
}
\value{
a list of data.table objects
}
\description{
Convert a list of fragment GenomicRanges to data.table
}
