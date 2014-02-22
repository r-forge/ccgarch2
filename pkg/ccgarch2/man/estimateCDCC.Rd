\name{estimateCDCC}
\alias{estimateCDCC}
\alias{print.cdcc}

\title{Estimating corrected DCC-GARCH model}
\description{This function estimates a corrected Dynamic Conditional Correlation GARCH 
        model.}
\usage{
estimateCDCC(inia = NULL, iniA = NULL, iniB = NULL, ini.dcc = NULL, 
      data, model="diagonal", \dots)

\method{print}{cdcc}(x, \dots)

}

\arguments{
  \item{inia}{a vector of initial values in the GARCH part.}
  \item{iniA}{a square or diagonal matrix of initial values for the ARCH parameter matrix in the GARCH part.}
  \item{iniB}{a square or diagonal matrix of initial values for the GARCH parameter matrix in the GARCH part.}
  \item{ini.dcc}{a vector of length 2 initial parameters for the cDCC part.}
  \item{data}{a data frame or a matrix object containing the variables.}
  \item{model}{a chacacter string setting the GARCH part of the model.}
  \item{x}{{an object of class \code{"cdcc"}.}}
  \item{...}{optional arguments. Currently not in use.}
}



\details{

  This function estimates a corrected Dynamic Conditional Correlation (cDCC-) GARCH 
  model of Aielli (2013), which is a modification of the original DCC-GARCH model 
  of Engle (2002). 
  
  The extractor function \code{summary()} is available for a \code{"cdcc"} class 
  object displaying a table with estimates and inferencial statistics, 
  information criterion and some diagnostic results of the standardized 
  residuals. See \code{\link{summary.cdcc}} for details.

  Estimation of the model is carried out in two steps. 
  The optimization in both steps is implemented by \code{\link[Rsolnp]{solnp}} function in 
  \pkg{Rsolnp} package and its outcome is save \sQuote{as is} in a list component 
  \code{f.stage} and \code{s.stage}, respectively.
  
  \code{\link[Rsolnp]{solnp}} uses a sequencial quadratic programming (SQP) technique 
  to optimize the objective function. See the manual and the references therein for 
  details. During the optimization, positivity and stationarity 
  restrictions are imposed on the parameters. 
  
  \pkg{Rsonlp} has a function called \code{\link[Rsolnp]{gosolnp}} that generates 
  initial starting values randomly. A logical argument \code{random} controlls which 
  function to use in the optimization. The default is \code{random = FALSE} meaning that 
  \code{solnp} of \pkg{Rsolnp} is used. If \code{random = TRUE}, \code{gosolnp} is used. 
  It should be noted that optimization takes longer time and may not always be successful 
  when \code{gosolnp} is used.

}

\value{

This function returns an S3 class object \code{"cdcc"} that is a list with the following components.

\item{f.stage}{an output of the first stage optimization from \pkg{Rsolnp}.}
\item{s.stage}{an output of the second stage optimization from \pkg{Rsolnp}.}
\item{model}{a chacacter string setting the GARCH part of the model.}
\item{method}{the optimization method used.}
\item{initial}{a list with the initial parameter vectors/matrices. The included components are \code{a}, \code{A}, \code{B}, and \code{dcc.par}.}
\item{data}{the original data matrix. Returned as a \pkg{zoo} object.}
\item{CDCC}{a matrix of the conditional correlations. Each row corresponds 
  to \eqn{\mathsf{vec}(\bm{R}_{t})}{vec(R_t)}. Returned as a \pkg{zoo} object.}
\item{h}{a matrix of the conditional variances. Each row corresponds 
  to \eqn{\bm{h}_{t}}{h_t}. Returned as a \pkg{zoo} object.}
\item{z}{a matrix of the standardized residuals. Each row corresponds 
  to \eqn{\bm{z}_{t}}{z_t}. Returned as a \pkg{zoo} object.}
}


\references{

  Aielli, G.P. (2013), 
  \dQuote{Dynamic Conditional Correlation: On Properties and Estimation.}
  \emph{Journal of Business and Economic Statistics}
  \bold{31}, 282--299.

  Engle, R.F. (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalized Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.

}

\seealso{
\code{\link{summary.cdcc}}, 
\code{\link{simulateCDCC}}, 
\code{\link{plot.cdcc}}, 
\code{\link[stats]{optim}} 
}

\examples{
\dontrun{
ndim <- 3
nobs <- 3000

## Setting parameters in the DCC part
a <- 0.05
b <- 0.8

## Setting a correlation matrix
Q <- diag(0, ndim)
Q[lower.tri(Q)] <- c(0.8, 0.3, 0.1)
Q <- Q + t(Q)
diag(Q) <- 1

## setting parameters in the GARCH part
a0 <- c(0.05, 0.07, 0.1)
A <- matrix(c(0.06, 0.0, 0.0, 0.07, 0.08, 0.002, 0.003, 0.0, 0.06), 
    ndim, ndim)
B <- matrix(c(0.75, 0.09, 0.03, 0.001, 0.81, 0.003, 0.001, 0.008, 0.84), 
    ndim, ndim)

## Simulating data
sim.data <- simulateCDCC(a, b, Q, a0, diag(diag(A)), diag(diag(B)), nobs)

## Estimating a cDCC model
garch_cdcc <- estimateCDCC(inia=a0, iniA=diag(diag(A)), iniB=diag(diag(B)), 
    ini.dcc=c(a, b), data=sim.data$eps, model="diagonal")

## Summarizing the results
summary_garch_cdcc <- summary(garch_cdcc)

## Plotting items
plot(garch_cdcc, item="correlation")   # for dynamic conditional correlations
plot(garch_cdcc, item="volatility")    # for volatility (square root of conditional variance)
plot(garch_cdcc, item="std.residuals") # for satandardized residuals
plot(garch_cdcc, item="return")        # for return (original data)
}
}


\keyword{ccgarch2}
\keyword{CDCC}
\keyword{nonlinear}

