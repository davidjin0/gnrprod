
prodgnr <- function(formula, share, id, time, data, in.price = 1,
                    out.price = 1,  ...) {
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "share", "id", "time", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- quote(stats::model.frame)

  mf <- eval.parent(mf)

  # Theorem 2 assumes access to rho/P (relative price deflator). Is it
  # reasonable to assume that the user has these parameters? Should we add
  # an option so that the user can add them? For now I set them to 1.
  # rho = 1
  # P = 1
  # s = ((rho * data[, vars[4]]) / (P * data[, vars[1]]))[, 1]
  # s = log(s)

  if (missing(share)) {
    mf <- mf %>% rename(
      y = names(mf)[1],
      l = names(mf)[2],
      k = names(mf)[3],
      i = names(mf)[4],
      id = names(mf)[5],
      time = names(mf)[6]
    )

    mf$s = log((in.price * mf$i) / (out.price *  mf*y))
  } else {
    mf <- mf %>% rename(
      y = names(mf)[1],
      l = names(mf)[2],
      k = names(mf)[3],
      i = names(mf)[4],
      s = names(mf)[5],
      id = names(mf)[6],
      time = names(mf)[7]
    )
  }

  # How to make it so that function can identify which inputs are flexible?
  # maybe add an option: organize(formula, data, intermediates = i)

  mf <- mf %>% mutate(
    y = log(y),
    l = log(l),
    k = log(k),
    i = log(i),
    ll = l * l,
    kk = k * k,
    ii = i * i,
    lk = l * k,
    li = l * i,
    ki = k * i,
    lki = l * k * i
  )

  mf <- na.omit(mf)

  first_stage_df <- first_stage(mf)
  point_estimates <- second_stage(first_stage_df)
  return(point_estimates)
}















