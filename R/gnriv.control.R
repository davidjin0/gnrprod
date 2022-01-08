#' Control iterations in second stage of GNR
#' @description Allows the user to modify convergence parameters of Gauss Newton algorithm used in the \code{gnriv} function
#' 
#' @param method temp
#' @param trace temp
#' @param fnscale temp
#' @param parscale temp
#' @param ndeps temp
#' @param maxit  temp
#' @param abstol temp
#' @param reltol temp
#' @param alpha temp
#' @param beta temp
#' @param gamma temp
#' @param REPORT temp
#' @param type temp
#' @param lmm temp
#' @param factr temp
#' @param pgtol temp
#' @param tmax temp
#' @param temp temp
#' @return a list of called control parameters
#'
#' @export

gnriv.control <- function(method = "BFGS", trace, fnscale, parscale, ndeps,
                          maxit, abstol, reltol, alpha, beta, gamma, REPORT,
                          type, lmm, factr, pgtol, tmax, temp) {
  return(list(method, match.call())[-2])
}
