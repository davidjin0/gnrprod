#' Estimate the nonparametric model in two stages: the share regression and instrumental variable estimation
#'
#' @param output name (character) of variable of level gross output or a numeric vector/dataframe/matrix
#' @param fixed name (character or character vector) of variables of level fixed inputs or a numeric vector/dataframe/matrix
#' @param flex name (character) of variable of level flexible input or a numeric vector/dataframe/matrix
#' @param share name (character) of variable of level intermediate input's revenue share or a numeric vector/dataframe/matrix
#' @param in_price optional name (character) of variable of common flexible input price or a numeric vector/dataframe/matrix
#' @param out_price optional name (character) of variable of common output price or a numeric vector/dataframe/matrix
#' @param id name (character) of variable of firm id or a numeric vector/dataframe/matrix
#' @param time name (character) of variable of time or a numeric vector/dataframe/matrix
#' @param data dataframe containing all variables with names specified by arguments above (left empty if arguments above are vector/dataframe/matrix)
#' @param degree degree of share regression polynomial
#' @param markov_degree degree of Markovian process for persistent productivity
#' @param firstStageOnly boolean for whether return is only first stage results
#' @param fs_control list of convergence parameters for first stage
#' @param ss_control list of convergence parameters for second stage
#' @return a list of class "gnr" with five elements:
#' \code{avg_elasticity}: a named numeric vector of the average elasticities of all inputs
#'
#' \code{data}: a list (data frame) containing:\code{output}, \code{fixed}, \code{flex}, \code{id}, \code{time}, and \code{share} variables and individual estimated elasticities
#'
#' \code{first_stage}: a list containing six elements describing the share regression (first stage):
#' \itemize{
#'  \item{\code{coef}}{: a numeric vector of the coefficients}
#'  \item{\code{residuals}}{: a numeric vector of the residuals}
#'  \item{\code{SSR}}{: sum of squared residual}
#'  \item{\code{iterations}}{: number of iterations performed}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved}
#'  \item{\code{control}}{: a list containing information on algorithm convergence (Gauss-Newton)}
#' }
#'
#' \code{second_stage}: a list containing four elements describing the second stage:
#' \itemize{
#'  \item{\code{productivity}}{: a numeric vector of the total productivity of each observations}
#'  \item{\code{iterations}}{: number of iterations performed}
#'  \item{\code{convergence}}{: boolean indicating whether convergence was achieved}
#'  \item{\code{control}}{: a list containing information on algorithm convergence (Gauss-Newton)}
#' }
#' @export

gnrprod <- function(output, fixed, flex, share, in_price = NULL,
                    out_price = NULL, id, time, data, degree = 2,
                    markov_degree = 2, firstStageOnly = FALSE,
                    fs_control = list(), ss_control = list()) {

  cl <- match.call()
  if (is.character(output)) {
    y <- as.matrix(log(data[, output]))
    output_name <- output
  } else {
    y <- as.matrix(log(output))
    output_name <- colnames(output)
    if (is.null(output_name)) {
      output_name <- deparse(substitute(output))
    }
  }

  if (is.character(fixed)) {
    fixed_input <- as.matrix(log(data[, fixed]))
    fixed_name <- fixed
  } else {
    fixed_input <- as.matrix(log(fixed))
    fixed_name <- colnames(fixed)
    if (is.null(fixed_name)) {
      fixed_name <- paste0(deparse(substitute(fixed)), 1L:ncol(fixed))
    }
  }

  if (is.character(flex)) {
    flex_input <- as.matrix(log(data[, flex]))
    flex_name <- flex
  } else {
    flex_input <- as.matrix(log(flex))
    flex_name <- colnames(flex)
    if (is.null(flex_name)) {
      flex_name <- deparse(substitute(flex))
    }
  }

  if (is.character(id)) {
    idm <- as.matrix(data[, id])
    id_name <- id
  } else {
    idm <- as.matrix(id)
    id_name <- colnames(id)
    if (is.null(id_name)) {
      id_name <- deparse(substitute(id))
    }
  }

  if (is.character(time)) {
    timem <- as.matrix(data[, time])
    time_name <- time
  } else {
    timem <- as.matrix(timem)
    time_name <- colnames(time)
    if (is.null(time_name)) {
      time_name <- deparse(substitute(time))
    }
  }

  if (!missing(share) && is.character(share)) {
    sharem <- as.matrix(log(data[, share]))
    share_name <- share
  } else if (!missing(share) && !is.character(share)) {
    sharem <- as.matrix(log(share))
    share_name <- colnames(share)
    if (is.null(share_name)) {
      share_name <- deparse(substitute(share))
    }
  } else if (!missing(in_price) && !missing(out_price)) {
    if (is.character(in_price)) {
      in_pricem <- as.matrix(data[, in_price])
    } else {
      in_pricem <- as.matrix(in_price)
    }

    if (is.character(out_price)) {
      out_pricem <- as.matrix(data[, out_price])
    } else {
      out_pricem <- as.matrix(out_price)
    }

    sharem = as.matrix(log((in_price * flex_input) / (out_price * y)))
    share_name = "share"
  } else {
    stop("must specify either share or intermediate-input price and output price")
  }

  complete_obs <- stats::complete.cases(cbind(y, fixed_input, flex_input, idm,
                                              timem, sharem))
  y <- y[complete_obs]
  fixed_input <- data.frame(fixed_input[complete_obs, ])
  flex_input <- data.frame(flex_input[complete_obs])
  idm <- idm[complete_obs]
  timem <- timem[complete_obs]
  sharem <- sharem[complete_obs]

  all_input = as.matrix(cbind(fixed_input, flex_input))
  poly_input = as.matrix(stats::poly(all_input, degree = degree, raw = TRUE))

  input_degrees <- sapply(colnames(poly_input), FUN = function(x) {
    a <- data.frame(base::strsplit(x, split = "[.]")[[1]])
  })
  input_degrees <- do.call(cbind, input_degrees)

  inputs = list("poly_input" = poly_input,
                "input_degrees" = input_degrees,
                "flex_in" = flex_input)

  first_stage_list <- first_stage(output = y,
                                input = inputs,
                                share = sharem,
                                id = idm,
                                time = timem,
                                fs_control = fs_control,
                                degree = degree)

  fs_results <- first_stage_list[[1]]
  fs_return <- first_stage_list[[2]]

  if (firstStageOnly) {
    return(fs_return)
  }

  second_stage_list <- second_stage(first_stage = fs_results, degree = degree,
                                    markov_degree = markov_degree,
                                    control = ss_control)

  ss_return <- second_stage_list[[2]]
  pred_elas <- second_stage_list[[1]]
  flex_elas <- fs_results$flex_in_elasticity

  elas <- data.frame(cbind(pred_elas, flex_elas))
  input_names <- c(fixed_name, flex_name)
  input_names <- paste(input_names, "elasticity", sep = "_")
  colnames(elas) <- input_names

  return_df <- data.frame(cbind(y, fixed_input, flex_input, idm, timem, sharem))
  colnames(return_df) <- c(output_name, fixed_name, flex_name, id_name,
                           time_name, share_name)

  point_estimates <- list("estimates" = elas, "data" = return_df)

  return_average_elas <- apply(elas, MARGIN = 2, FUN = mean)
  names(return_average_elas) <- paste(input_names, "avg", sep = "_")
  gnr_out <- list(return_average_elas, point_estimates, fs_return, ss_return, cl)
  names(gnr_out) <- c("avg_elasticity", "point_estimates", "first_stage",
                          "second_stage", "call")
  class(gnr_out) <- "gnr"
  return(gnr_out)
}



