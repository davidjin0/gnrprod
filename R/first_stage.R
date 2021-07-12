library(tidyverse)

rm(list = ls())
first_stage <- function(formula, data) {
  vars = all.vars(formula)
  s <- data[, vars[1]]
  l <- log(data[, vars[2]])
  k <- log(data[, vars[3]])
  i <- log(data[, vars[4]])

  df <- data.frame(s, l, k, i)

  df <- df %>% rename(
    s = si,
    l = L,
    k = K,
    i = RI
  )

  df <- df %>% mutate(
    ll = l * l,
    kk = k * k,
    ii = i * i,
    lk = l * k,
    li = l * i,
    ki = k * i,
    lki = l * k * i
  )

  df <- na.omit(df)

  start_reg = stats::lm(s ~ l + k + i + ll + kk + ii + lk + li + ki + lki,
                        data = df)

  constant <- predict(start_reg, df) - coef(start_reg)[1]
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

  estimates <- gauss_newton_reg(start, df, 1)
  return(estimates)
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

  new_start <- start + (inv(t(X) %*% X) %*% t(X) %*% initial_errors)
  new_pred <- evaluate(new_start, data)
  new_errors <- cbind(data$s - log(new_pred))
  new_SSR <- t(new_errors) %*% new_errors

  while (is.na(new_SSR) || new_SSR > initial_SSR) {
    alpha <- alpha / 2
    new_start <- start +
      alpha * (inv(t(X) %*% X) %*% t(X) %*% initial_errors)
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

data = read_dta('data_col.dta')
firststage <- first_stage(si ~ L + K + RI, data)







