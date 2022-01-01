
second_stage <- function(first_stage, degree, markov_degree, control) {
  in_df <- data.frame(rbind(first_stage$input, first_stage$input_degree))
  pred <- in_df[, in_df[nrow(in_df), ] == 0]
  pred <- pred[-c((nrow(pred) -
                     (nrow(first_stage$input_degree) - 1)):(nrow(pred))), ]
  flex <- in_df[, in_df[nrow(in_df), ] != 0]
  flex <- flex[-c((nrow(flex) -
                     (nrow(first_stage$input_degree) - 1)):(nrow(flex))), ]
  id <- data.frame(first_stage$id)
  time <- data.frame(first_stage$time)
  big_Y <- data.frame(first_stage$big_Y)
  errors <- first_stage$errors

  constant_reg <- stats::lm(as.matrix(big_Y) ~ as.matrix(pred))
  coefficients <- as.matrix(stats::coef(constant_reg)[2:(ncol(pred) + 1)])

  ss_df <- data.frame(cbind(id, time, pred, big_Y))
  ss_df <- lag_gnr(ss_df)

  pred_base <- as.matrix(ss_df[[1]])
  big_Y_base <- pred_base[, ncol(pred_base)]
  pred_base <- pred_base[, -c(ncol(pred_base))]

  pred_lag <- as.matrix(ss_df[[2]])
  big_Y_lag <- pred_lag[, ncol(pred_lag)]
  pred_lag <- pred_lag[, -c(ncol(pred_lag))]

  ctrl <- do.call("second_stage_control", control)
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }

  constant_gmm <- pracma::gaussNewton(x0 = coefficients, data = pred_base,
                              Ffun = constant_moments,
                              big_Y_base = big_Y_base,
                              big_Y_lag = big_Y_lag,
                              lag_data = pred_lag,
                              markov_degree = markov_degree,
                              maxiter = ctrl$maxit,
                              tol = ctrl$reltol)

  convergence <- ifelse(constant_gmm$niter >= ctrl$maxit, FALSE, TRUE)
  C_coef <- constant_gmm$xs

  input_degree <- first_stage$input_degree
  all_input <- first_stage$input
  gamma <- first_stage$gamma[-1, ]
  gamma_denom <- first_stage$gamma_denom[-1, ]

  elasticities <- lapply(1:(nrow(input_degree) - ncol(first_stage$flex_base)),
                         FUN = function(i) {
    new_in_deg <- input_degree
    new_in_deg[i, ] <- ifelse(new_in_deg[i, ] > 0,
                              new_in_deg[i, ] - 1,
                              new_in_deg[i, ])

    new_C_deg <- new_in_deg[, new_in_deg[nrow(input_degree), ] == 0]

    new_in_deg[nrow(new_in_deg), ] <- new_in_deg[nrow(new_in_deg), ] + 1

    in_match <- apply(new_in_deg, MARGIN = 2, FUN = match_gnr, degree_vec =
                        input_degree)
    C_match <- apply(new_C_deg, MARGIN = 2, FUN = match_gnr, degree_vec =
                       input_degree)

    deriv_input <- all_input[, in_match]
    deriv_input[is.na(deriv_input)] <- 0

    deriv_C <- all_input[, C_match]
    deriv_C[is.na(deriv_C)] <- 1

    C <- deriv_C %*%
      t(t(input_degree[i, input_degree[nrow(input_degree), ] == 0]) * C_coef)

    elas <- deriv_input %*% t(t(input_degree[i, ]) * (gamma / gamma_denom)) + C
  })

  logomega <- big_Y - (as.matrix(pred) %*% C_coef)
  omega <- exp(logomega)
  productivity <- as.matrix(exp(logomega + errors))

  ss_return <- list(productivity, markov_degree, constant_gmm$niter, ctrl$maxit,
                    ctrl$reltol, convergence)
  names(ss_return) <- c("productivity", "markov_degree", "iterations", "maxit",
                        "reltol", "convergence")

  average_elas <- lapply(elasticities, FUN = mean)
  elasticities <- do.call(cbind, elasticities)
  return_list <- list(elasticities, ss_return)
  return(return_list)
}

constant_moments <- function(C_kl, data, big_Y_base, big_Y_lag, lag_data,
                             markov_degree) {
  w = big_Y_base - (data %*% C_kl)
  w_1 = big_Y_lag - (lag_data %*% C_kl)

  poly = sapply(2:markov_degree, FUN = function(i) {
    `^`(w_1, i)
  })

  markov = cbind(w_1, poly)

  reg <- stats::lm(w ~ markov)
  csi <- w - reg$fitted.values

  moments <- apply(data, MARGIN = 2, FUN = function(i) {
    sum(i * csi) / length(i)
  })

  return(moments)
}

match_gnr <- function(i, degree_vec) {
  match_test <- apply(i == degree_vec, MARGIN = 2, FUN = prod)
  col_index <- ifelse(length(which(match_test == 1)) != 0,
                      which(match_test == 1), NA)
}

second_stage_control <- function(maxit = 100, reltol = 1e-8) {
  if (maxit <= 0 || !is.numeric(maxit)) {
    stop("maximum iteration count must be integer > 0")
  }

  if (reltol <= 0 || !is.numeric(reltol)) {
    stop("relative tolerance must be numeric > 0")
  }

  list(maxit = maxit, reltol = reltol)
}



