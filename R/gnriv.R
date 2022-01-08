#' Estimate fixed input elasticity and total productivity: Gandhi, Navarro, Rivers (GNR) lag instruments
#' @description The \code{gnriv} function implements the second stage of the GNR production function estimation routine, nonparametrically identifying the fixed input elasticities of the production function and total productivity. This function is called in the main wrapper function \code{gnrprod}. It currently supports only one flexible input.
#'
#' For details, see Gandhi, Navarro, and Rivers (2020).
#'
#' @param object object of class \code{gnrflex}.
#' @param degree degree of Markov process for persistent productivity. Defaults to 2.
#' @param control an optional list of convergence settings. See \code{gnriv.control} for listing.
#' @param ... additional optional arguments passed to optim.
#' @return a list of class "gnriv" containing six elements:
#'
#' \code{pred_elas}: a numeric matrix of estimated elasticities of fixed inputs for each observation.
#'
#' \code{productivity}: a numeric vector of estimated total productivity.
#'
#' \code{degree}: degree of Markov process.
#'
#' \code{iterations}: number of iterations performed.
#'
#' \code{convergence}: boolean indicating whether convergence was achieved.
#'
#' \code{control}: the list of convergence control parameters. See \code{gnriv.control}.
#' 
#' @usage gnriv(object, degree = 2, control, ...)
#' 
#' @export


gnriv <- function(object, degree = 2, control, ...) {
  if (attr(object, "class") != "gnrflex") {
    stop("object must be of class gnrflex")
  }

  in_df <- data.frame(rbind(object$arg$input, object$arg$input_degree))
  pred <- in_df[, in_df[nrow(in_df), ] == 0]
  pred <- pred[-c((nrow(pred) -
                     (nrow(object$arg$input_degree) - 1)):(nrow(pred))), ]
  flex <- in_df[, in_df[nrow(in_df), ] != 0]
  flex <- flex[-c((nrow(flex) -
                     (nrow(object$arg$input_degree) - 1)):(nrow(flex))), ]

  id <- object$arg$id
  time <- object$arg$time
  big_Y <- object$arg$big_Y
  errors <- object$elas$residuals

  constant_reg <- stats::lm(big_Y ~ as.matrix(pred))
  coefficients <- stats::coef(constant_reg)[2:(ncol(pred) + 1)]

  ss_df <- data.frame(cbind(id, time, pred, big_Y))
  ss_df <- lag_gnr(ss_df)

  pred_base <- as.matrix(ss_df[[1]])
  big_Y_base <- pred_base[, ncol(pred_base)]
  pred_base <- pred_base[, -c(ncol(pred_base))]

  pred_lag <- as.matrix(ss_df[[2]])
  big_Y_lag <- pred_lag[, ncol(pred_lag)]
  pred_lag <- pred_lag[, -c(ncol(pred_lag))]

  ctrl <- gnriv.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }
  method = ctrl[[1]]
  ctrl[[1]] = NULL

  constant_gmm <- stats::optim(par = coefficients, fn = constant_moments,
                               data = pred_base, big_Y_base = big_Y_base,
                               big_Y_lag = big_Y_lag, lag_data = pred_lag,
                               degree = degree, method = method,
                               control = ctrl, ...)
  
  opt_ctrl <- list(trace = 0, fnscale = 1,
                   parscale = rep.int(1, length(coefficients)),
                   ndeps = rep.int(1e-3, length(coefficients)),
                   maxit = 100L, abstol = -Inf,
                   reltol = sqrt(.Machine$double.eps),
                   alpha = 1.0, beta = 0.5, gamma = 2.0,
                   REPORT = 10,
                   type = 1,
                   lmm = 5, factr = 1e7, pgtol = 0,
                   tmax = 10, temp = 10.0)
  
  opt_ctrl[names(ctrl)] <- ctrl[names(ctrl)]
  C_coef <- constant_gmm$par

  input_degree <- object$arg$input_degree
  all_input <- object$arg$input

  elasticities <- lapply(1:(nrow(input_degree) - 1), FUN = function(i) {
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

    elas <- deriv_input %*% t(t(input_degree[i, ]) * (object$arg$D_coef)) + C
  })

  logomega <- big_Y - (as.matrix(pred) %*% C_coef)
  omega <- exp(logomega)
  productivity <- as.matrix(exp(logomega + errors))

  elasticities <- do.call(cbind, elasticities)
  ss_return <- list("pred_elas" = elasticities,
                    "productivity" = productivity,
                    "degree" = degree,
                    "method" = method,
                    "opt_info" = constant_gmm,
                    "control" = opt_ctrl)
  class(ss_return) <- "gnriv"
  return(ss_return)
}

constant_moments <- function(C_kl, data, big_Y_base, big_Y_lag, lag_data,
                             degree) {
  w <- big_Y_base - (data %*% C_kl)
  w_1 <- big_Y_lag - (lag_data %*% C_kl)

  poly <- sapply(2:degree, FUN = function(i) {
    `^`(w_1, i)
  })

  markov <- cbind(w_1, poly)

  reg <- stats::lm(w ~ markov)
  csi <- w - reg$fitted.values

  moments <- apply(data, MARGIN = 2, FUN = function(i) {
    sum(i * csi) / length(i)
  })
  
  obj <- t(moments) %*% moments
  
  return(obj)
}

match_gnr <- function(i, degree_vec) {
  match_test <- apply(i == degree_vec, MARGIN = 2, FUN = prod)
  col_index <- ifelse(length(which(match_test == 1)) != 0,
                      which(match_test == 1), NA)
}


