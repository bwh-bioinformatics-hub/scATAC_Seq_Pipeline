#' Compare two matrices, and find the column in the second matrix with maximum correlation to each column in the first matrix
#'
#' @param query_mat A matrix object to use as a query.
#' @param target_mat A matrix object to use as a target.
#' @param method A character object specifying the correlation method to use. Passed to base R cor(). Default is "pearson".
#' @param n_threads A numeric object specifying the number of threads to use. Default is 1.
#'
#' @return a data.frame with stats for the top-2 best correlated columns in target for each column of query.
#' @export
max_column_cor <- function (query_mat,
                            target_mat,
                            method = "pearson",
                            n_threads = 1) {
  if(n_threads > 1) {
    out_list <- parallel::mclapply(1:ncol(query_mat),
                                   function(i) {
                                     query_vals <- query_mat[, i]
                                     cor_vals <- apply(target_mat, 2, function(x) {
                                       cor(query_vals, x, method = method)
                                     })

                                     cor_vals <- cor_vals[order(cor_vals, decreasing = TRUE)]

                                     data.frame(query = colnames(query_mat)[i],
                                                max_cor = cor_vals[1],
                                                max_target = names(cor_vals)[1],
                                                second_cor = cor_vals[2],
                                                second_target = names(cor_vals)[2],
                                                max_separation = cor_vals[1] - cor_vals[2])
                                   },
                                   mc.cores = n_threads)
  } else {

    out_list <- lapply(1:ncol(query_mat),
                       function(i) {
                         query_vals <- query_mat[, i]
                         cor_vals <- apply(target_mat, 2, function(x) {
                           cor(query_vals, x, method = method)
                         })

                         cor_vals <- cor_vals[order(cor_vals, decreasing = TRUE)]

                         data.frame(query = colnames(query_mat)[i],
                                    max_cor = cor_vals[1],
                                    max_target = names(cor_vals)[1],
                                    second_cor = cor_vals[2],
                                    second_target = names(cor_vals)[2],
                                    max_separation = cor_vals[1] - cor_vals[2])
                       })
  }

  out_df <- do.call("rbind",
                    out_list)
}
