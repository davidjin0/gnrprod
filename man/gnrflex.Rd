% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gnrflex.R
\name{gnrflex}
\alias{gnrflex}
\title{Estimate flexible input elasticity: Gandhi, Navarro, Rivers (GNR) share regression; first stage}
\usage{
gnrflex(output, fixed, flex, share, id, time, data, control)
}
\arguments{
\item{output}{name (character) of variable of log gross output in data or a numeric vector.}

\item{fixed}{name (character or character vector) of variables of log fixed inputs in data or a numeric matrix.}

\item{flex}{name (character) of variable of log flexible input in data or a numeric vector.}

\item{share}{name (character) of variable of log intermediate input's revenue share in data or a numeric vector.}

\item{id}{name (character) of variable of firm id in data or a numeric vector.}

\item{time}{name (character) of variable of time in data or a numeric vector.}

\item{data}{\code{\link[base]{data.frame}} containing all variables with names specified by arguments above (left empty if arguments above are vector/matrix rather than strings).}

\item{control}{an optional list of convergence settings. See \code{\link[gnrprod]{gnrflex.control}} for listing.}
}
\value{
a list of class 'gnrflex' containing three elements:

\code{elas}: a list containing six elements describing the share regression:
\itemize{
\item{\code{flex_elas}}{: a numeric vector of the estimated flexible input elasticity for each observation.}
\item{\code{coefficients}}{: a numeric vector of the coefficients of the estimator scaled by a constant. See Gandhi, Navarro, and Rivers (2020, p. 2994, equation (21)).}
\item{\code{residuals}}{: a numeric vector of the residuals.}
\item{\code{SSR}}{: sum of squared residuals.}
\item{\code{iterations}}{: number of iterations performed.}
\item{\code{convergence}}{: boolean indicating whether convergence was achieved.}
}

\code{arg}: a list containing eight elements to be passed to the second stage function \code{\link[gnrprod]{gnriv}}:
\itemize{
\item{\code{input}}{: a numeric matrix (S3: \code{\link[stats]{poly}}) of the polynomial expansion of all inputs.}
\item{\code{input_degree}}{: a numeric matrix corresponding to \code{input} denoting each vector's degree.}
\item{\code{all_input}}{: a numeric matrix of the inputs without polynomial expansion.}
\item{\code{big_Y}}{: a numeric vector of persistent productivity minus the constant of integration. See Gandhi, Navarro, and Rivers (2020, p. 2991, equation (16)).}
\item{\code{D_coef}}{: a numeric vector equaling \code{coef} divided by an estimate of the constant.}
\item{\code{id}}{: a numeric vector of the firm ids.}
\item{\code{time}}{: a numeric vector of time.}
\item{\code{degree}}{: the degree of the share regression.}
\item{\code{fixed_names}}{: the names of fixed inputs. To be used in the second stage.}
}

\code{control}: the list of convergence control parameters. See \code{\link[gnrprod]{gnrflex.control}} for available parameters.
}
\description{
The \code{gnrflex} function implements the first stage (share
regression) of the GNR production function estimation routine,
nonparametrically identifying the flexible input elasticity of the
production function. This function is called within the main wrapper
function \code{\link[gnrprod]{gnrprod}}. If the production-related inputs
are characters, a \code{\link[base]{data.frame}} must be specified under
\code{data}. Alternatively, matrices/vectors may be directly specified
without specifying \code{data}. \code{gnrprod} currently supports only one
flexible input. The parameters are optimized using the Gauss-Newton
algorithm. \code{gnrflex} currently supports only one flexible input.

For details, see Gandhi, Navarro, and Rivers (2020).
}
\examples{
require(gnrprod)
data <- colombian
industry_311_flex <- gnrflex(output = "RGO", fixed = c("L", "K"),
                             flex = "RI", share = "share", id = "id",
                             time = "year", data = data,
                             control = list(degree = 2, maxit = 200))
}
\references{
Gandhi, Amit, Salvador Navarro, and David Rivers. 2020. "On the Identification of Gross Output Production Functions." \emph{Journal of Political Economy}, 128(8): 2973-3016. \doi{10.1086/707736}.

Davidson, Russell, James G. MacKinnon. 1993. "The Gauss-Newton Regression." In \emph{Estimation and Inference in Econometrics}, 176-207. New York: Oxford University Press.
}
