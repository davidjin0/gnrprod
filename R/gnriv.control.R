#' Control parameters in \code{gnriv}
#' @description Allows the user to modify convergence parameters of Gauss Newton algorithm used in the \code{\link[gnrprod]{gnriv}} function
#' 
#' @param degree degree of Markov process for persistent productivity. Defaults to 3.
#' @param method the method of optimization passed to \code{\link[stats]{optim}}. Defaults to "BFGS." See \code{\link[stats]{optim}} under 'Details' for listing of available methods.
#' @param ... additional optional control parameters passed to \code{\link[stats]{optim}}. See \code{\link[stats]{optim}} for available parameters.
#' @return a list containing \code{degree} and \code{method} and any additional parameters in \code{...}. 
#' 
#' @usage gnriv.control(degree = 3, method = "BFGS", ...)
#' 
#' @export

gnriv.control <- function(degree = 3, method = "BFGS", ...) {
  return(list(degree = degree, method = method, ...))
}



