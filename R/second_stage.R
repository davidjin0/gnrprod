
second_stage <- function(data) {
  data <- data %>% group_by(id) %>%
    arrange(time, .by_group = TRUE) %>%
    mutate(big_Y_1 = dplyr::lag(big_Y, n = 1, default = NA),
    l_1 = dplyr::lag(l, n = 1, default = NA),
    k_1 = dplyr::lag(k, n = 1, default = NA),
    ll_1 = dplyr::lag(ll, n = 1, default = NA),
    kk_1 = dplyr::lag(kk, n = 1, default = NA),
    lk_1 = dplyr::lag(lk, n = 1, default = NA),
    big_Y_2 = dplyr::lag(big_Y_1, n = 1, default = NA),
    l_2 = dplyr::lag(l_1, n = 1, default = NA),
    k_2 = dplyr::lag(k_1, n = 1, default = NA),
    ll_2 = dplyr::lag(ll_1, n = 1, default = NA),
    kk_2 = dplyr::lag(kk_1, n = 1, default = NA),
    lk_2 = dplyr::lag(lk_1, n = 1, default = NA))

  constant_reg <- lm(big_Y ~ l + k + ll + kk + lk, data = data)

  coefficients <- (as.matrix(coef(constant_reg)))[2:6, ]

  al <- as.numeric(coefficients[1])
  ak <- as.numeric(coefficients[2])
  al2 <- as.numeric(coefficients[3])
  ak2 <- as.numeric(coefficients[4])
  alk <- as.numeric(coefficients[5])

  df_gmm <- complete_obs(data, cols = c("big_Y", "l", "k", "big_Y_1", "l_1", "k_1"))

  df_gmm <- df_gmm %>% select(
    big_Y, l, k, ll, kk, lk, big_Y_1, l_1, k_1, ll_1, kk_1, lk_1, id, time
  )

  constant_gmm <- gaussNewton(x0 = c(al, ak, al2, ak2, alk), data = df_gmm,
                              Ffun = constant_moments)

  C_coef = constant_gmm$xs

  data <- data %>% mutate(
    logomega = big_Y - (C_coef[1] * l) - (C_coef[2] * k) - (C_coef[3] * ll) -
      (C_coef[4] * kk) - (C_coef[5] * lk),
    omega = exp(logomega),
    prod = exp(logomega + errors),

    l_elas = (gl * i + 2 * gll * li + glk * ki) + (gli * ii / 2) + (glki * ki * i/2)
    + C_coef[1] + (2 * C_coef[3] * l) + (C_coef[5] * k),
    k_elas = (gk * i + 2 * gkk * ki + glk * li) + (gki * ii / 2) + (glki * li * i/2)
    + C_coef[2] + (2 * C_coef[4] * k) + (C_coef[5] * l)
  )

  mean_lelas <- mean(data$l_elas, na.rm = T)
  mean_kelas <- mean(data$k_elas, na.rm = T)
  mean_ielas <- mean(data$i_elas, na.rm = T)
  mean_sum_elas = mean(data$l_elas + data$k_elas + data$i_elas, na.rm = T)
  mean_l_k_elas = mean_kelas / mean_lelas

  point_estimates <- c(mean_lelas, mean_kelas, mean_ielas, mean_sum_elas,
                       mean_l_k_elas)
  return(point_estimates)
}

complete_obs <- function(data, cols) {
  complete_vectors <- complete.cases(data[, cols])
  return(data[complete_vectors, ])
}

constant_moments <- function(C_kl, data) {
  data <- data %>% mutate(
    w = big_Y - {C_kl[1] * l} - {C_kl[2] * k} - {C_kl[3] * ll} - {C_kl[4] * kk} - {C_kl[5] * lk},
    w_1 = big_Y_1 - {C_kl[1] * l_1} - {C_kl[2] * k_1} - {C_kl[3] * ll_1} - {C_kl[4] * kk_1} - {C_kl[5] * lk_1},
    w2_1 = w_1 * w_1,
    w3_1 = w_1 * w_1 * w_1
  )

  reg <- stats::lm(w ~ w_1 + w2_1 + w3_1, data)
  data$csi <- data$w - predict(reg, data)

  m1 <- sum(data$l * data$csi) / nrow(data)
  m2 <- sum(data$k * data$csi) / nrow(data)
  m3 <- sum(data$ll * data$csi) / nrow(data)
  m4 <- sum(data$kk * data$csi) / nrow(data)
  m5 <- sum(data$lk * data$csi) / nrow(data)

  M <- c(m1, m2, m3, m4, m5)
  return(M)
}






