coef.gnr <- function(object) object$avg_elasticity

print.gnr <- function(object, digits = max(3L, getOption("digits") - 3L)) {
  cat("Call:\n")
  cat("  output: ", deparse(object$call$output), "\n", sep = "")
  cat("  fixed inputs: ", deparse(object$call$fixed), "\n", sep = " ")
  cat("  flexible inputs: ", deparse(object$call$flex), "\n", sep = " ")
  if (!is.null(object$call$data)) {
    cat("  data: ", deparse(object$call$data), "\n", sep = "")
  }
  cat("\n")
  cat("Average input elasticities:\n")
  print(object$avg_elasticity, digits = digits)
  cat("\n")
  cat("First stage sum of squared residuals: ", object$first_stage$SSR)
  invisible(object)
}

extract.gnr <-function(object) {
  cbind(object$data, object$second_stage$productivity)
}