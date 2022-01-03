#' Control iterations in second stage of GNR
#' @description Allows the user to modify convergence parameters of Gauss Newton algorithm used in the \code{gnriv} function
#'
#' @param maxit Maximum number of iterations. Defaults to 100.
#' @param reltol Relative convergence tolerance. Defaults to 1e-8.
#' @return a list containing two elements: \code{maxit} and \code{reltol}
#'
#' @export

gnriv.control <- function(maxit = 100, reltol = 1e-8) {
  if (maxit <= 0 || !is.numeric(maxit)) {
    stop("maximum iteration count must be integer > 0")
  }

  if (reltol <= 0 || !is.numeric(reltol)) {
    stop("relative tolerance must be numeric > 0")
  }

  list(maxit = maxit, reltol = reltol)
}
