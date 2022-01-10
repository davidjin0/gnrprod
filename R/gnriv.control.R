#' Control iterations in second stage of GNR
#' @description Allows the user to modify convergence parameters of Gauss Newton algorithm used in the \code{gnriv} function
#' 
#' @param degree degree of Markov process for persistent productivity. Defaults to 3.
#' @param method the method of optimization passed to optim. Defaults to "BFGS." See optim for listing of available methods.
#' @param ... additional optional control parameters passed to optim. See optim for available parameters.
#' @return a list containing `degree` and `method` and any additional parameters in `...`. 
#'
#' @export

gnriv.control <- function(degree = 3, method = "BFGS", ...) {
  
  return(list(degree = degree, method = method, ...))
}



