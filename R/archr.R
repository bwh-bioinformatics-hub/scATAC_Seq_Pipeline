#' Get ArchR matrices in dgCMatrix format
#'
#' @param proj an ArchRProject object
#' @param target The name of a matrix in the ArchRProject. Options are: "GeneScoreMatrix", "PeakMatrix", "TileMatrix"
#' @param ... Additional parameters passed to ArchR::getMatrixFromProject()
#'
#' @return a dgCMatrix object
#' @export
#'
get_archr_dgCMatrix <- function(proj,
                                target = "GeneScoreMatrix",
                                ...) {

  if(target == "GeneScoreMatrix") {

    se <- ArchR::getMatrixFromProject(
      proj,
      useMatrix = "GeneScoreMatrix",
      ...)

    mat <- se@assays@data@listData$GeneScoreMatrix

    rownames(mat) <- se@elementMetadata@listData$name

  } else if(target == "PeakMatrix") {
    se <- ArchR::getMatrixFromProject(
      proj,
      useMatrix = "PeakMatrix",
      ...)

    mat <- se@assays@data@listData$PeakMatrix

    rownames(mat) <- paste(
      seqnames(se@rowRanges),
      start(se@rowRanges),
      end(se@rowRanges),
      sep = "_")

  } else if(target == "TileMatrix") {
    se <- ArchR::getMatrixFromProject(
      proj,
      useMatrix = "TileMatrix",
      binarize = TRUE,
      ...)

    mat <- se@assays@data@listData$TileMatrix

    rownames(mat) <- paste(
      seqnames(se@rowRanges),
      start(se@rowRanges),
      end(se@rowRanges),
      sep = "_")
  }

  mat

}

#' Get Peak-to-Gene links from an ArchRProject
#'
#' @param proj an ArchRProject object
#' @param fdr_cut cutoff for p2g FDR. Default is 0.05
#' @param cor_cut cutoff for p2g correlation. Default is 0.3.
#' @param cor_type one of "pos" for positive correlations, "neg" for negative correlations, or "any" for either direction. Default is "any".
#'
#' @return a data.frame of Peak-to-Gene links with FDR and correlation scores.
#' @export
#'
get_archr_p2g_df <- function(proj,
                                fdr_cut = 0.05,
                                cor_cut = 0.3,
                                cor_type = "any") {

  peak_meta <- proj@peakSet@metadata

  p2g_df <- as.data.frame(peak_meta$Peak2GeneLinks@listData)
  p2g_df <- p2g_df[!is.na(p2g_df$Correlation),]
  p2g_df <- p2g_df[p2g_df$FDR < fdr_cut,]

  if(cor_type == "pos") {
    p2g_df <- p2g_df[p2g_df$Correlation > abs(cor_cut),]
  } else if(cor_type == "neg") {
    p2g_df <- p2g_df[p2g_df$Correlation < -1 * abs(cor_cut),]
  } else if(cor_type == "any") {
    p2g_df <- p2g_df[abs(p2g_df$Correlation) > abs(cor_cut),]
  }

  p2g_df
}

#' Pull the peakSet object from an ArchRProject
#'
#' @param proj an ArchRProject object
#'
#' @return a peakSet object
#' @export
#'
get_archr_peakSet <- function(proj) {
  peak_meta <- proj@peakSet@metadata
  peak_meta$Peak2GeneLinks@metadata$peakSet
}

#' Pull the geneSet object from an ArchRProject
#'
#' @param proj an ArchRProject object
#'
#' @return a geneSet object
#' @export
#'
get_archr_geneSet <- function(proj) {
  peak_meta <- proj@peakSet@metadata
  peak_meta$Peak2GeneLinks@metadata$geneSet
}

#' Get Peak-to-Gene correlations from an ArchRProject as a sparse matrix
#'
#' @param proj an ArchRProject object
#' @param fdr_cut cutoff for p2g FDR. Default is 0.05
#' @param cor_cut cutoff for p2g correlation. Default is 0.3.
#' @param cor_type one of "pos" for positive correlations, "neg" for negative correlations, or "any" for either direction. Default is "any".
#'
#' @return a Matrix::dgCMatrix of Peak-to-Gene links with correlation scores as values
#' @export
#'
get_archr_p2g_mat <- function(proj,
                              fdr_cut = 0.05,
                              cor_cut = 0.3,
                              cor_type = "any") {

  peak_set <- get_archr_peakSet(proj)
  gene_set <- get_archr_geneSet(proj)

  p2g_df <- get_archr_p2g_df(proj,
                             fdr_cut = fdr_cut,
                             cor_cut = cor_cut,
                             cor_type = cor_type)

  Matrix::sparseMatrix(
    i = p2g_df$idxATAC,
    j = p2g_df$idxRNA,
    x = p2g_df$Correlation,
    index1 = TRUE,
    dims = c(length(peak_set), length(gene_set)),
    dimnames = list(paste0("peak_",1:length(peak_set)),
                    gene_set@elementMetadata@listData$name)
  )
}

#' Get a Peak Annotation sparse matrix from an ArchRProject
#'
#' @param proj an ArchRProject object
#' @param peakAnnotation The annotation to retrieve, e.g. "EncodeTFBS" or "Motif"
#'
#' @return a sparse matrix (usually logical) indicating which peaks match which annotations
#' @export
#'
get_archr_peakAnno_mat <- function(proj, peakAnnotation = "EncodeTFBS") {
  anno_file <- proj@peakAnnotation@listData[[peakAnnotation]]$Matches
  anno_res <- readRDS(anno_file)

  assays(anno_res)@listData$matches
}

#' Link peakAnno to genes via Peak-to-Gene linkage
#'
#' This is matrix multiplication of the peakAnnotation matrix with the Peak-to-Gene matrix
#'
#' @param anno_mat a peakAnnotation matrix, as retrieved with get_archr_peakAnno_mat()
#' @param p2g_mat a Peak-to-Gene matrix, as retrieved with get_archr_p2g_mat()
#'
#' @return a Matrix::dgCMatrix of peakAnno-to-Gene links
#' @export
#'
link_peakAnno_to_gene <- function(anno_mat, p2g_mat) {
  Matrix::t(anno_mat) %*% p2g_mat
}

#' Convert a sparse linkage matrix to an edge data.frame for graph analysis
#'
#' @param link_mat a matrix with sources as columns, targets as rows, and weights as values.
#'
#' @return a data.frame with columns "from", "to", and "weight".
#' @export
#'
link_mat_to_df <- function(link_mat) {
  data.frame(
    from = rep(colnames(link_mat), diff(link_mat@p)),
    to = rownames(link_mat)[link_mat@i + 1],
    weight = link_mat@x
  )
}

#' Convert a sparse linkage matrix to a graph object
#'
#' @param link_mat a matrix with sources as columns, targets as rows, and weights as values.
#'
#' @return an igraph object
#' @export
#'
link_mat_to_graph <- function(link_mat) {
  tidygraph::as_tbl_graph(
    link_mat_to_df(link_mat)
  )
}


chr_fragments_datatable <- function(h5_con, chr_name) {
  rg <- h5read(h5_con, paste0("/Fragments/",chr_name,"/Ranges"))

  data.table(
    chr = chr_name,
    start = rg[,1],
    end = rg[,1] + rg[,2],
    bc = rep(h5read(h5_con, paste0("/Fragments/",chr_name,"/RGValues")),
             h5read(h5_con, paste0("/Fragments/",chr_name,"/RGLengths"))),
    n_umi = 1
  )
}

extract_arrow_fragments <- function(arrow_file) {
  h5_con <- H5Fopen(arrow_file)
  contents <- h5ls(h5_con)
  chrs <- contents$name[contents$group == "/Fragments"]

  frags <- chr_fragments_datatable(h5_con, chrs[1])
  for(i in 2:length(chrs)) {
    frags <- rbind(frags, chr_fragments_datatable(h5_con, chrs[i]))
  }

  setorderv(frags, c("chr","start"))

  frags
}

read_arrow_cell_meta <- function(arrow_file) {
  h5_con <- H5Fopen(arrow_file)
  contents <- h5ls(h5_con)

  meta_names <- contents[contents$group == "/Metadata",]
  meta_names <- meta_names$name[meta_names$dim != 1]

  meta_list <- lapply(
    meta_names,
    function(meta_name) {
      vals <- h5read(h5_con, paste0("/Metadata/",meta_name))
      as.vector(vals)
    }
  )

  df <- as.data.frame(meta_list)
  names(df) <- meta_names

  df
}
