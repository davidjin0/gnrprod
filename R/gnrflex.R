#' Estimate flexible input elasticity: Gandhi, Navarro, Rivers (GNR) share regression
#' @description The \code{gnrflex} function implements the first stage (share regression) of the GNR production function estimation routine, nonparametrically identifying the flexible input elasticity of the production function. This function is called in the main wrapper function \code{gnrprod}. It currently supports only one flexible input.
#'
#' For details, see Gandhi, Navarro, and Rivers (2020).
#'
#' @param output a numeric vector of level gross output
#' @param fixed a numeric matrix of level fixed inputs
#' @param flex a numeric vector of level flexible input
#' @param share a numeric vector of the level intermediate input's revenue share
#' @param id a numeric vector of firm ids
#' @param time a numeric vector of time
#' @param degree degree of share regression polynomial
#' @param control an optional list of convergence settings. See \code{gnrflex.control} for listing.
#' @return a list of class "gnrflex" containing three elements:
#'
#' \code{elas}: a list containing six elements describing the share regression:
#'
#' \itemize{
#'  \item{\code{flex_in_elas}}{: a numeric vector of the estimated flexible input elasticity for each observation}
#'  \item{\code{coef}}{: a numeric vector of the coefficients of the estimator scaled by a constant (equation (21))}
#'  \item{\code{residuals}}{: a numeric vector of the residuals}
#'  \item{\code{SSR}}{: sum of squared residuals}
#'  \item{\code{iterations}}{: number of iterations performed}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved}
#' }
#'
#' \code{arg}: a list containing seven elements to be passed to the second stage function \code{gnriv}
#' \itemize{
#'  \item{\code{input}}{: a numeric matrix (S3 'poly') of the polynomial expansion of all inputs}
#'  \item{\code{input_degree}}{: a numeric matrix corresponding to \code{input} denoting each vector's degree}
#'  \item{\code{big_Y}}{: a numeric vector of persistent productivity minus the constant of integration (equation (16) in Gandhi, Navarro, and Rivers (2020))}
#'  \item{\code{D_coef}}{: a numeric vector equalling \code{coef} divided by an estimate of the constant}
#'  \item{\code{id}}{: a numeric vector of the firm ids}
#'  \item{\code{time}}{: a numeric vector of time}
#'  \item{\code{degree}}{: the degree of the share regression}
#' }
#'
#' \code{control}: the list of convergence control parameters. See \code{gnrflex.control}
#'
#' @export

gnrflex <- function(output, fixed, flex, share, id, time, degree = 2,
                    control = gnrflex.control()) {
  
  complete_obs <- stats::complete.cases(cbind(output, fixed, flex, id, time,
                                              share))
  
  if (length(complete_obs) != length(output)) {
    stop("'output', 'fixed', 'flex', 'share', 'id', and 'time' cannot contain missing values")
  }
  
  all_input <- cbind(fixed, flex)
  poly_input <- stats::poly(all_input, degree = degree, raw = TRUE)
  input_degrees <- sapply(colnames(poly_input), FUN = function(x) {
    a <- base::strsplit(x, split = "[.]")[[1]]
  })
  input_degrees <- apply(input_degrees, 2, as.numeric)

  gamma_denom <- rbind(1, as.matrix(input_degrees[nrow(input_degrees), ] + 1))
  start_reg <- stats::lm(share ~ poly_input)

  constant <- start_reg$fitted.values - stats::coef(start_reg)[1]
  constant <- -min(constant, na.rm = TRUE) + 0.1

  start <- c(constant, (coef(start_reg)[-1]))
  share_reg <- gauss_newton_reg(start = start, data = poly_input,
                                share = share, control = control)

  coef <- share_reg[[1]]
  i_elas <- log(pred_fs(coef, poly_input))
  errors <- i_elas - share
  mean_exp_err <- mean(exp(errors))
  i_elas <- exp(i_elas - log(mean_exp_err))
  gamma <- as.matrix(coef / mean_exp_err)
  flex_gamma <- gamma / gamma_denom

  integ_G_I <- pred_fs(flex_gamma, poly_input)
  integ_G_I <- integ_G_I * flex
  big_Y <- as.matrix(output - errors - integ_G_I)

  fs_elas <- list("flex_elas" = i_elas,
                  "coef" = coef,
                  "residuals" = errors,
                  "SSR" = share_reg$SSR,
                  "iterations" = share_reg$iterations,
                  "convergence" = share_reg$convergence)

  fs_arg <- list("input" = poly_input,
                 "input_degree" = input_degrees,
                 "big_Y" = big_Y,
                 "D_coef" = gamma[-1, ] / gamma_denom[-1, ],
                 "id" = id,
                 "time" = time,
                 "degree" = degree)

  fs_return <- list("elas" = fs_elas, "arg" = fs_arg, "control" = control)
  class(fs_return) <- "gnrflex"
  return(fs_return)
}

gauss_newton_reg <- function(start, data, share, control) {
  ctrl <- gnrflex.control()
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }

  iter <- 0
  call_start <- start

  inputs_1 <- data.frame(rep(1, nrow(data)), data)
  names(inputs_1)[1] <- "constant"

  while (iter < ctrl$maxit) {
    initial_pred <- pred_fs(call_start, data)

    X <- as.matrix(inputs_1 / initial_pred)

    initial_errors <- as.matrix(cbind(share - log(initial_pred)))
    initial_SSR <- t(initial_errors) %*% initial_errors

    new_start <- call_start + (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
    new_pred <- pred_fs(new_start, data)
    suppressWarnings(new_errors <- cbind(share - log(new_pred)))
    new_SSR <- t(new_errors) %*% new_errors

    initial_step <- ctrl$initial_step
    min_factor <- ctrl$min_factor

    while ((is.na(new_SSR) || new_SSR > initial_SSR)
           & initial_step >= min_factor) {
      initial_step <- initial_step / 2
      new_start <- call_start +
        initial_step * (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
      new_pred <- pred_fs(new_start, data)
      suppressWarnings(new_errors <- cbind(share - log(new_pred)))
      new_SSR <- t(new_errors) %*% new_errors
    }

    conv_bool <- TRUE
    for (i in 1:length(new_start)) {
      if (initial_step * abs(new_start[i]) > ctrl$reltol *
          (abs(new_start[i]) + 1e-3)) {
        conv_bool <- FALSE
        break
      }
    }

    iter <- iter + 1
    if (conv_bool) {
      return_list <- list(new_start, iter, new_SSR, conv_bool, control)
      names(return_list) <- c("share_reg_coef", "iterations", "SSR",
                              "convergence")
      return(return_list)
    }
    call_start <- new_start
  }
  warning("share regression failed to converge")
  return_list <- list(new_start, iter, new_SSR, conv_bool)
  names(return_list) <- c("share_reg_coef", "iterations", "SSR", "convergence")
  return(return_list)
}


pred_fs <- function(start, data) {
  matrix <- as.matrix(cbind(rep(1, nrow(data)), data))
  new_m <- matrix %*% start
  return(new_m)
}



