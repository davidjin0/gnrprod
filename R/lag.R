lag_gnr <- function(unsorted) {
  unique_ids <- unique(unsorted[, 1])
  lag_indic <- 3:ncol(unsorted)

  lagged_list <- lapply(unique_ids, FUN = function(i) {
    i_df <- unsorted[unsorted[, 1] == i, ]
    i_df <- i_df[order(i_df[, 2]), ]

    lag_1_df <- data.frame(matrix(NA, nrow = 1, ncol = length(lag_indic)))
    colnames(lag_1_df)[1:ncol(lag_1_df)] = colnames(i_df)[lag_indic]
    lag_1_df <- rbind(lag_1_df, i_df[, lag_indic])
    lag_1_df <- lag_1_df[-nrow(lag_1_df), ]

    i_df <- cbind(i_df, lag_1_df)
  })

  lagged_df <- do.call(rbind, lagged_list)
  lagged_df <- stats::na.omit(lagged_df)

  lag_base_list <- list(lagged_df[, 3:ncol(unsorted)],
                        lagged_df[, (ncol(unsorted) + 1):ncol(lagged_df)])

  return(lag_base_list)
}









