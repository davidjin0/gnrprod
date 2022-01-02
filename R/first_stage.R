
first_stage <- function(output, input, share, id, time, fs_control,
                        degree) {
  all_input <- as.matrix(input$poly_input)
  input_degrees <- data.frame(input$input_degrees)
  input_degrees <- sapply(input_degrees, as.numeric)

  gamma_denom <- rbind(1, matrix(input_degrees[nrow(input_degrees), ] + 1))
  start_reg <- stats::lm(share ~ all_input)

  constant <- start_reg$fitted.values - stats::coef(start_reg)[1]
  constant <- -min(constant, na.rm = T) + 0.1

  start <- c(constant, (coef(start_reg)[-1]))

  share_reg <- gauss_newton_reg(start = start, data = all_input, share = share,
                                control = fs_control)

  coef = share_reg[[1]]
  i_elas = log(evaluate_gnr(coef, all_input))

  errors = i_elas - share
  mean_exp_err = mean(exp(errors))
  i_elas = exp(i_elas - log(mean_exp_err))

  gamma <- as.matrix(coef / mean_exp_err)
  flex_gamma <- gamma / gamma_denom

  integ_G_I = evaluate_gnr(flex_gamma, all_input)
  integ_G_I = integ_G_I * input$flex_in
  big_Y = output - errors - integ_G_I

  first_stage_results <- list("input" = all_input,
                              "flex_base" = input$flex_in,
                              "input_degree" = input_degrees,
                              "id" = id,
                              "time" = time,
                              "big_Y" = big_Y,
                              "errors" = errors,
                              "flex_in_elasticity" = i_elas,
                              "gamma" = gamma,
                              "gamma_denom" = gamma_denom)

  first_stage_return <- list("coef" = coef,
                             "residuals" = errors,
                             "SSR" = share_reg$SSR,
                             "iterations" = share_reg$iterations,
                             "convergence" = share_reg$convergence)
  return(list(first_stage_results, first_stage_return))
}



gauss_newton_reg <- function(start, data, share, control) {
  ctrl <- do.call("first_stage_control", control)
  if (!missing(control)) {
    control <- as.list(control)
    ctrl[names(control)] <- control
  }

  iter <- 0
  call_start <- start

  inputs_1 = data.frame(rep(1, nrow(data)), data)
  names(inputs_1)[1] <- "constant"

  while (iter < ctrl$maxit) {
    initial_pred <- evaluate_gnr(call_start, data)

    X <- as.matrix(inputs_1 / initial_pred)

    initial_errors <- as.matrix(cbind(share - log(initial_pred)))
    initial_SSR <- t(initial_errors) %*% initial_errors

    new_start <- call_start + (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
    new_pred <- evaluate_gnr(new_start, data)
    suppressWarnings(new_errors <- cbind(share - log(new_pred)))
    new_SSR <- t(new_errors) %*% new_errors

    initial_step = ctrl$initial_step
    min_factor = ctrl$min_factor

    while ((is.na(new_SSR) || new_SSR > initial_SSR)
           & initial_step >= min_factor) {
      initial_step <- initial_step / 2
      new_start <- call_start +
        initial_step * (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
      new_pred <- evaluate_gnr(new_start, data)
      suppressWarnings(new_errors <- cbind(share - log(new_pred)))
      new_SSR = t(new_errors) %*% new_errors
    }

    conv_bool <- TRUE
    for (i in 1:length(new_start)) {
      if (initial_step * abs(new_start[i]) > ctrl$reltol *
          (abs(new_start[i]) + 1e-3)) {
        conv_bool <- FALSE
        break
      }
    }

    iter = iter + 1
    if (conv_bool) {
      return_list <- list(new_start, iter, new_SSR, conv_bool)
      names(return_list) <- c("share_reg_coef", "iterations", "SSR",
                              "convergence")
      return(return_list)
    }
    call_start = new_start
  }
  warning("share regression failed to converge")
  return_list <- list(new_start, iter, new_SSR, conv_bool)
  names(return_list) <- c("share_reg_coef", "iterations", "SSR", "convergence")
  return(return_list)
}


first_stage_control <- function(maxit = 100, reltol = 1e-5, initial_step = 100,
                                min_factor = 1e-5) {
  if (maxit <= 0 || !is.numeric(maxit)) {
    stop("maximum iteration count must be integer > 0")
  }

  if (reltol <= 0 || !is.numeric(reltol)) {
    stop("relative tolerance must be numeric > 0")
  }

  if (initial_step <= 0 || !is.numeric(initial_step)) {
    stop("initial_step must be numeric > 0")
  }

  if (min_factor < 0 || !is.numeric(min_factor)) {
    stop("min_factor must be numeric >= 0")
  }

  list(maxit = maxit, reltol = reltol, initial_step = initial_step,
       min_factor = min_factor)
}


evaluate_gnr <- function(start, data) {
  matrix = as.matrix(cbind(rep(1, nrow(data)), data))
  new_m = matrix %*% start
  return(new_m)
}





