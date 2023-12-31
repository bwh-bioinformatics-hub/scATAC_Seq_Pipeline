% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{window_index_to_bed}
\alias{window_index_to_bed}
\title{Convert a vector of window indexes to a BED-like table of chromosome positions for each window.}
\usage{
window_index_to_bed(window_indexes, chrom_sizes, window_size = 5000)
}
\arguments{
\item{window_indexes}{a numeric vector of window indexes}

\item{chrom_sizes}{A data.frame of chromosome sizes, e.g. as generated by read_chrom_sizes()}

\item{window_size}{A numeric value specifying the size of windows to use as features. Default is 5000.}
}
\value{
a data.table of chromosome positions, with columns chr, start, and end
}
\description{
Convert a vector of window indexes to a BED-like table of chromosome positions for each window.
}
