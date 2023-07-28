#' Vectorized Chi-squared tests for differential gene/accessibility detection
#'
#' This function uses vectors of the number of samples in two sets that have detection of
#' a set of genes and the total number of cells in each set to compute Chi-quared tests with 1 DOF for
#' differential detection.
#'
#' @param x an integer vector with the number of cells in group \emph{x} with detection of each gene.
#' @param x.total an integer value with the total number of cells in group \emph{x}.
#' @param y an integer vector with the number of cells in group \emph{y} with detection of each gene.
#' @param y.total an integer value with the total number of cells in group \emph{y}.
#'
#' @return a data.frame with the following result for each gene:
#' \itemize{
#' \item{stats: The value of the chi-squared test statistic}
#' \item{pval: The p-value as reported by pchisq}
#' \item{logFC: The log2(fold change) in detection frequency between samples (x / y)}
#' \item{diff: The difference in proportions between the samples (x - y)}
#' }
#'
#' @export
vec_chisq_test <- function(x,
                           x.total,
                           y,
                           y.total) {

  total <- x.total + y.total
  present <- x + y
  absent <- total - x - y

  o <- cbind(x,
             x.total - x,
             y,
             y.total - y)

  e <- cbind(present * x.total,
             absent * x.total,
             present * y.total,
             absent * y.total)

  e <- e / as.vector(total)

  stat <- rowSums(pmax(0, abs(o - e) - 0.5) ^ 2 / e)

  results <- data.frame(stats = stat,
                        pval = pchisq(stat, 1, lower.tail = FALSE),
                        logFC = log2( (x * y.total) / (y * x.total)),
                        diff = x / x.total - y / y.total)
  results
}

#' Chi-squared tests for rows of a sparse matrix
#'
#' @param mat A sparse matrix with observations as columns and features as rows.
#' @param fg_cols Column names for foreground columns
#' @param bg_cols Column names for background columns. If NULL (default) uses all non-fg_cols.
#' @param rows Row names to use for performing tests. If NULL (default), uses all rows
#' @param sort_on A character object specifying which column to use for sorting results. If NULL, will not sort. Default is "padj".
#'
#' @return a data.frame of Chi-squared test results.
#' @export
#'
sparse_mat_chisq <- function(mat,
                             fg_cols,
                             bg_cols = NULL,
                             rows = NULL,
                             sort_on = "padj") {

  assertthat::assert_that(sum(fg_cols %in% colnames(mat)) == length(fg_cols))

  if(!is.null(rows)) {
    mat <- mat[rows,]
  }

  if(is.null(bg_cols)) {
    bg_cols <- setdiff(colnames(mat), fg_cols)
  } else {
    assertthat::assert_that(sum(bg_cols %in% colnames(mat)) == length(fg_cols))
  }

  mat@x <- rep(1, length(mat@x))

  fg_counts <- rowSums(mat[,fg_cols])
  bg_counts <- rowSums(mat[,bg_cols])

  chisq_res <- vec_chisq_test(x = fg_counts,
                              x.total = length(fg_cols),
                              y = bg_counts,
                              y.total = length(bg_cols))

  chisq_res <- cbind(feature = rownames(mat),
                     chisq_res)
  chisq_res$padj <- p.adjust(chisq_res$pval,
                             method = "BH")

  if(!is.null(sort_on)) {
    chisq_res <- chisq_res[order(chisq_res[[sort_on]]),]
  }

  chisq_res
}

