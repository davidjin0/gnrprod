
first_stage <- function(data) {
  start_reg = stats::lm(s ~ l + k + i + ll + kk + ii + lk + li + ki + lki,
                        data = data)

  constant <- predict(start_reg, data) - coef(start_reg)[1]
  constant <- -min(constant, na.rm = T) + 0.1

  start <- c(
    constant,
    (coef(start_reg)["l"]),
    (coef(start_reg)["k"]),
    (coef(start_reg)["i"]),
    (coef(start_reg)["ll"]),
    (coef(start_reg)["lk"]),
    (coef(start_reg)["li"]),
    (coef(start_reg)["kk"]),
    (coef(start_reg)["ki"]),
    (coef(start_reg)["ii"]),
    (coef(start_reg)["lki"])
  )

  share_reg <- gauss_newton_reg(start, data, 1)

  data$i_elas = log(evaluate(share_reg, data))

  data <- data %>% dplyr::mutate(
    errors = i_elas - s,
    mean_exp_err = mean(exp(errors)),
    i_elas = exp(i_elas - log(mean_exp_err))
  )

  data <- data %>% mutate(
    g0 = share_reg[1] / mean_exp_err,
    gl = share_reg[2] / mean_exp_err,
    gk = share_reg[3] / mean_exp_err,
    gi = share_reg[4] / mean_exp_err,
    gll = share_reg[5] / mean_exp_err,
    glk = share_reg[6] / mean_exp_err,
    gli = share_reg[7] / mean_exp_err,
    gkk = share_reg[8] / mean_exp_err,
    gki = share_reg[9] / mean_exp_err,
    gii = share_reg[10] / mean_exp_err,
    glki = share_reg[11] / mean_exp_err
  )

  data <- data %>% mutate(
    integ_G_I = (g0 + (gl * l) + (gk * k) + (gll * ll) + (gkk * kk) +
                   (glk * lk) + (gi * i + gli * li + gki * ki + glki * lki) / 2 +
                   (gii * ii) / 3) * i,
    big_Y = y - errors - integ_G_I
  )

  return(data)
}



gauss_newton_reg <- function(start, data, alpha) {
  initial_pred <- evaluate(start, data)

  X <- cbind(
    1 / initial_pred,
    data$l / initial_pred,
    data$k / initial_pred,
    data$i / initial_pred,
    data$ll / initial_pred,
    data$lk / initial_pred,
    data$li / initial_pred,
    data$kk / initial_pred,
    data$ki / initial_pred,
    data$ii / initial_pred,
    data$lki / initial_pred
  )

  initial_errors <- cbind(data$s - log(initial_pred))
  initial_SSR <- t(initial_errors) %*% initial_errors

  new_start <- start + (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
  new_pred <- evaluate(new_start, data)
  new_errors <- cbind(data$s - log(new_pred))
  new_SSR <- t(new_errors) %*% new_errors

  while (is.na(new_SSR) || new_SSR > initial_SSR) {
    alpha <- alpha / 2
    new_start <- start +
      alpha * (solve(t(X) %*% X) %*% t(X) %*% initial_errors)
    new_pred <- evaluate(new_start, data)
    new_errors <- cbind(data$s - log(new_pred))
    new_SSR = t(new_errors) %*% new_errors
  }

  conv_bool <- TRUE
  for (i in 1:11) {
    if (alpha * abs(new_start[i]) > 1e-5 * (abs(new_start[i]) + 1e-3)) {
      conv_bool <- FALSE
      break
    }
  }

  if (conv_bool) {
    return(new_start)
  }

  return(gauss_newton_reg(new_start, data, 1))
}

evaluate <- function(start, data) {
  value <- start[1] + (start[2] * data$l) + (start[3] * data$k) +
    (start[4] * data$i) + (start[5] * data$ll) + (start[6] * data$lk) +
    (start[7] * data$li) + (start[8] * data$kk) + (start[9] * data$ki) +
    (start[10] * data$ii) + (start[11] * data$lki)

  return(value)
}






