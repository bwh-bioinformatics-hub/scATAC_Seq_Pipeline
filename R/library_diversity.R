
#' Project library diversity using preseqR
#'
#' @param freq_mat a matrix containing the histogram of counts
#' @param total_metrics a data.frame of sequencing metrics
#' @param step_size a numeric value specifying how large of a step to use for projection. Default is 1e7
#' @param max_val a numeric value specifying the maximum extent of projection. Default is 5e9
#' @param times a numeric value specifying how many rounds of bootstrapping to use. Passed to preseqR.rSAC.bootstrap()
#'
#' @return a data.frame of projection results.
#' @export
diversity_projection <- function(freq_mat,
                                 total_metrics,
                                 step_size = 1e7,
                                 max_val = 5e9,
                                 times = 100) {

  est <- preseqR::preseqR.rSAC.bootstrap(n = freq_mat,
                                         times = times)

  umis_per_step <- floor(step_size * total_metrics$total_umis / total_metrics$total_reads)
  frac_per_step <- umis_per_step / total_metrics$total_umis

  max_frac <- max_val / step_size

  steps <- seq(from = 0,
               to = frac_per_step * max_frac,
               by = frac_per_step)

  proj <- suppressWarnings(
    data.frame(n_raw_reads = total_metrics$total_reads * steps,
               n_mapped_reads = total_metrics$total_counts * steps,
               expected_umis = est$f(steps),
               ci_95_low = est$lb(steps),
               ci_95_high = est$ub(steps))
  )

  proj

}
