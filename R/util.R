#' Binarize/convert all values in a sparse matrix to 1
#'
#' @param mat a sparse matrix object from the Matrix package, e.g. dgCMatrix.
#'
#' @return a sparse matrix object
#' @export
binarize <- function(mat) {
  mat@x <- rep(1, length(mat@x))
  mat
}


#' Load a stored reference dataset:
#'
#' @param name one of:
#' \itemize{
#'   \item{"altius"}{ENCODE/Altius index reference regions from Mueleman, et al. 2020 (ENCFF503GCK)}
#'   \item{"gene_bodies"}{Gene body regions from Ensembl (hg38 = v93; hg19 = v87)}
#'   \item{"great"}{Gene-adjacent regions computed a la GREAT}
#'   \item{"peaks"}{PBMC peaks identified in Lareau, et al. 2019 (GSE123577)}
#'   \item{"tss"}{TSS +/- 2kb regions from Ensembl (hg38 = v93; hg19 = v87)}
#' }
#' @param genome either "hg38" or "hg19".
#' @param format either "gr" for GenomicRanges or "dt" for data.table
#'
#' @return output based on the format parameter
#' @export
#'
load_reference <- function(name = c("altius","gene_bodies","great","peaks","tss"),
                           genome = c("hg38","hg19"),
                           format = c("gr","dt")) {

  name <- match.arg(name)
  genome <- match.arg(genome)
  format <- match.arg(format)

  if(format == "gr") {
    suffix <- "gr.rds"
  } else if(format == "dt") {
    suffix <- "bed.gz"
  }

  filename <- paste(genome, name, suffix, sep = "_")
  ref_path <- system.file("reference", package = "ATACSeqPipeline")

  if(format == "gr") {
    readRDS(file.path(ref_path, filename))
  } else if(format == "dt") {
    data.table::fread(file.path(ref_path, filename))
  }

}
