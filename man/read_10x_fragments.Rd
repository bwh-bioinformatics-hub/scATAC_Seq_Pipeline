% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/io.R
\name{read_10x_fragments}
\alias{read_10x_fragments}
\title{Read 10x scATAC-seq fragment data}
\usage{
read_10x_fragments(
  fragments_tsv,
  singlecell_csv,
  min_reads = 1000,
  remove_chrM = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{fragments_tsv}{A fragments.tsv.gz file generated by cellranger-atac count in the outs/ directory}

\item{singlecell_csv}{A singlecell.csv file generated by cellranger-atac count in the outs/ directory for filtering. If NULL, will read all fragments.}

\item{min_reads}{A numeric value indicating the minimum number of reads required to retain a cell barcode. Default is 1000.}

\item{remove_chrM}{A logical value indicating whether or not to remove all fragments mapping to the mitochondrial genome. Default is TRUE}

\item{verbose}{A logical value indicating whether or not to display messages. Default is TRUE.}
}
\value{
a list of GenomicRanges objects, one per cell barcode
}
\description{
Read 10x scATAC-seq fragment data
}
