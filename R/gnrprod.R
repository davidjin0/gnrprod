#' Estimate production functions and productivity: Gandhi, Navarro, and Rivers (2020)
#'
#' @description The \code{gnrprod} function is the front end of the \code{gnrprod} package. It estimates production functions and productivity in two stages: \code{gnrflex} (estimate flexible input elasticity) and \code{gnriv} (estimate fixed input elasticities and productivity). It currently supports only one flexible input.
#'
#' @param output name (character) of variable of level gross output or a numeric vector
#' @param fixed name (character or character vector) of variables of level fixed inputs or a numeric matrix
#' @param flex name (character) of variable of level flexible input or a numeric vector
#' @param share name (character) of variable of level intermediate input's revenue share or a numeric vector
#' @param in_price optional (required if \code{share} is not specified) name (character) of variable of common flexible input price or a numeric vector
#' @param out_price optional (required if \code{share} is not specified) name (character) of variable of common output price or a numeric vector
#' @param id name (character) of variable of firm id or a numeric vector
#' @param time name (character) of variable of time or a numeric vector
#' @param data dataframe containing all variables with names specified by arguments above (left empty if arguments above are vector/matrix)
#' @param degree degree of share regression polynomial
#' @param markov_degree degree of Markov process for persistent productivity
#' @param fs_control an optional list of convergence settings of the first stage. See \code{gnrflex.control} for listing.
#' @param ss_control an optional list of convergence settings of the second stage. See \code{gnriv.control} for listing.
#' @return a list of class "gnr" with five elements:
#' \code{avg_elasticity}: a named numeric vector of the average elasticities of all inputs
#'
#' \code{data}: a list (dataframe) containing:\code{output}, \code{fixed}, \code{flex}, \code{id}, \code{time}, and \code{share} variables and estimated elasticities for each observation
#'
#' \code{first_stage}: a list containing six elements describing the share regression (first stage):
#' \itemize{
#'  \item{\code{coef}}{: a numeric vector of the coefficients of the estimator scaled by a constant (equation (21))}
#'  \item{\code{residuals}}{: a numeric vector of the residuals}
#'  \item{\code{SSR}}{: sum of squared residual}
#'  \item{\code{iterations}}{: number of iterations performed}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved}
#'  \item{\code{control}}{: list of convergence control parameters (see \code{gnrflex.control})}
#' }
#'
#' \code{second_stage}: a list containing four elements describing the second stage:
#' \itemize{
#'  \item{\code{productivity}}{: a numeric vector of the estimated total productivity}
#'  \item{\code{iterations}}{: number of iterations performed}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved}
#'  \item{\code{control}}{: list of convergence control parameters (see \code{gnriv.control})}
#' }
#'
#' @usage gnrprod(output, fixed, flex, share, in_price = NULL,
#'                out_price = NULL, id, time, data, degree = 2,
#'                markov_degree = 2, fs_control = gnrflex.control(),
#'                ss_control = gnriv.control())
#' 
#' @examples
#' require(gnrprod)
#' data <- colombian
#' industry_311 <- gnrprod(output = "RGO", fixed = c("L", "K"), flex = "RI",
#'                         share = "share", id = "id", time = "year",
#'                         data = data, degree = 3, markov_degree = 3)
#' @export

gnrprod <- function(output, fixed, flex, share, in_price = NULL,
                    out_price = NULL, id, time, data, degree = 2,
                    markov_degree = 2, B = NULL, fs_control = gnrflex.control(),
                    ss_control = gnriv.control()) {

  cl <- match.call()
  
  output <- get_matrix(output)
  fixed <- get_matrix(fixed)
  flex <- get_matrix(flex)
  id <- get_matrix(id)
  time <- get_matrix(time)
  
  if (!missing(share)) {
    share <- get_matrix(share)
  } else if (!missing(in_price) && !missing(out_price)) {
    share <- (in_price + flex) - (out_price + output)
    in_price <- get_matrix(in_price)
    out_price <- get_matrix(out_price)
    colnames(share) <- "share1"
  } else {
    stop("must specify either share or both intermediate-input price and output price")
  }

  complete_obs <- stats::complete.cases(cbind(output, fixed, flex, share, id,
                                              time))
  output <- output[complete_obs, , drop = FALSE]
  fixed <- fixed[complete_obs, , drop = FALSE]
  flex <- flex[complete_obs, , drop = FALSE]
  id <- id[complete_obs, , drop = FALSE]
  time <- time[complete_obs, , drop = FALSE]
  share <- share[complete_obs, , drop = FALSE]
  
  if (!is.null(in_price) && !is.null(out_price)) {
    in_price <- in_price[complete_obs, , drop = FALSE]
    out_price <- out_price[complete_obs, , drop = FALSE]
    mf <- cbind(output, fixed, flex, in_price, out_price, share,
                           id, time)
  } else {
    mf <- cbind(output, fixed, flex, share, id, time)
  }
  
  mf <- mf[order(mf[, ncol(mf) - 1], mf[, ncol(mf)]), , drop = FALSE]
  
  
  gnr_flex <- gnrflex(output = output, fixed = fixed, flex = flex,
                      share = share, id = id, time = time, degree = degree,
                      control = fs_control)
  
  
  
  time_start <- Sys.time()
  gnr_iv <- gnriv(object = gnr_flex, degree = markov_degree,
                  control = ss_control)
  time_end <- Sys.time() - time_start
  return(gnr_iv)
  
  boot_sd = NULL
  if (!missing(B) && B > 1) {
    id_unique <- unique(id)
    boot_elas <- lapply(1:B, FUN = function(i) {
      boot_ids <- sample(id_unique, length(id_unique), replace = TRUE)
      boot_df <- do.call(rbind, lapply(boot_ids, function (x) {
        mf[mf[, colnames(id)] == x,]
      }))
      
      boot_output <- boot_df[, 1]
      boot_fixed <- boot_df[, 2:(ncol(fixed) + 1)]
      boot_flex <- boot_df[, (ncol(fixed) + 2):(ncol(fixed) + 1 + ncol(flex))]
      boot_share <- boot_df[, ncol(boot_df) - 2]
      boot_id <- boot_df[, ncol(boot_df) - 1]
      boot_time <- boot_df[, ncol(boot_df)]
      
      boot_fs <- gnrflex(output = boot_output, fixed = boot_fixed,
                         flex = boot_flex, share = boot_share, id = boot_id,
                         time = boot_time, degree = degree, control = fs_control)
      
      flex_elas <- mean(boot_fs$elas$flex_elas)
      
      boot_ss <- gnriv(object = boot_fs, degree = markov_degree,
                       control = ss_control)
      
      fixed_elas <- apply(boot_ss$elas, 2, mean)
      
      return(c(fixed_elas, flex_elas))
    })
    boot_est <- do.call(cbind, boot_elas)
    boot_sd <- apply(boot_est, 1, sd)
  }
  
  pred_elas <- gnr_iv$elas
  flex_elas <- gnr_flex$elas$flex_elas
  elas <- data.frame(cbind(pred_elas, flex_elas))
  input_names <- c(colnames(fixed), colnames(flex))
  input_names <- paste(input_names, "elasticity", sep = "_")
  colnames(elas) <- input_names
  mf <- cbind(mf, elas, gnr_iv$productivity, gnr_flex$elas$residuals)
  colnames(mf)[(ncol(mf) - 1):ncol(mf)] = c("productivity", "flex_resid")

  fs_return <- list("coefficients" = gnr_flex$elas$coef,
                    "SSR" = gnr_flex$elas$SSR,
                    "iterations" = gnr_flex$elas$iterations,
                    "convergence" = gnr_flex$elas$convergence,
                    "control" = gnr_flex$control)

  ss_return <- list("iterations" = gnr_iv$iterations,
                    "convergence" = gnr_iv$convergence,
                    "control" = gnr_iv$control)

  return_average_elas <- apply(elas, MARGIN = 2, FUN = mean)
  names(return_average_elas) <- paste(input_names, "avg", sep = "_")
  
  return_sd <- boot_sd
  if (!is.null(boot_sd)) {
    names(return_sd) <- paste(input_names, "se", sep = "_")
  }
  param = list("elas" = return_average_elas, "std_errors" = return_sd)
  
  gnr_out <- list(param, mf, fs_return, ss_return, cl)
  names(gnr_out) <- c("estimates", "data", "first_stage", "second_stage",
                      "call")
  class(gnr_out) <- "gnr"
  return(gnr_out)
}

get_matrix <- function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  
  if (is.character(x)) {
    col <- as.matrix(data[, x])
    colnames(col) <- x
  } else {
    col <- as.matrix(x)
    colnames(col) <- colnames(x, do.NULL = FALSE, prefix = paste(substitute(x)))
  }
  return(col)
}






