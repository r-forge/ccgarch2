\name{estimateCCC}
\alias{estimateCCC}
\alias{print.ccc}


\title{Estimating CCC-GARCH model}
\description{This function estimates a Constant Conditional Correlation GARCH 
        model.}
\usage{
estimateCCC(inia = NULL, iniA = NULL, iniB = NULL, 
    data, model="diagonal", \dots)

\method{print}{ccc}(x, \dots)
}

\arguments{
  \item{inia}{a vector of initial values in the GARCH part.}
  \item{iniA}{a square or diagonal matrix of initial values for the ARCH parameter matrix in the GARCH part.}
  \item{iniB}{a square or diagonal matrix of initial values for the GARCH parameter matrix in the GARCH part.}
  \item{data}{a data frame or a matrix object containing the variables.}
  \item{model}{a chacacter string setting the GARCH part of the model.}
  \item{x}{{an object of class \code{"dcc"}.}}
  \item{...}{optional arguments. Currently not in use.}
}

\details{
  This function estimates a Constant Conditional Correlation (CCC-) GARCH 
  model of Bollerslev (1990). 

  The extractor function \code{summary()} is available for a \code{"ccc"} class 
  object displaying a table of estimates and inferencial statistics, 
  information criterion and some diagnostic results of the standardized 
  residuals. See \code{\link{summary.ccc}} for details.

  Estimation of the model is carried out in a single step, that is, the parameters 
  in the GARCH part and the conditional correlations are simultaneously estimated. 
  The optimization is implemented by \code{\link[Rsolnp]{solnp}} function in 
  \pkg{Rsolnp} package and its outcome is save \sQuote{as is} in a list component 
  \code{results}.
  
  \code{\link[Rsolnp]{solnp}} uses a sequencial quadratic programming (SQP) technique 
  to optimize the objective function. See the manual and the references therein for 
  details. During the optimization, positivity and stationarity 
  restrictions are imposed on the parameters. 

}

\value{
This function returns an S3 class object \code{"ccc"} that is a list with the following components.

\item{results}{an output of the optimization from \code{optim}.}
\item{model}{a chacacter string setting the GARCH part of the model.}
\item{method}{the optimization method used.}
\item{initial}{a list with the initial parameter vectors/matrices. The included components are \code{a}, \code{A}, and \code{B}.}
\item{data}{the data matrix. Returned as a \pkg{zoo} object.}
\item{estimates}{a list with the following components: \code{a} for constant, 
  \code{A} for ARCH, \code{B} for GARCH parameters in the conditional variance equation, 
  and \code{R} for the constant conditional correlation matrix.}
\item{h}{a matrix of the conditional variances. Each row corresponds 
  to \eqn{\bm{h}_{t}}{h_t}. Returned as a \pkg{zoo} object.}
\item{z}{a matrix of the standardized residuals. Each row corresponds 
  to \eqn{\bm{z}_{t}}{z_t}. Returned as a \pkg{zoo} object.}
}


\references{
  Bollerslev, T. (1990),
  \dQuote{Modeling the Coherence in Short-Run Nominal Exchange 
          Rates: A Multivariate Generalized ARCH Approach},
  \emph{Review of Economics and Statistics},
  \bold{72}, 498--505.

  Engle, R.F. and K. Sheppard (2001),
  \dQuote{Theoretical and Empirical Properties of Dynamic
  Conditional Correlation Multivariate GARCH.}
  \emph{Stern Finance Working Paper Series}
  FIN-01-027 (Revised in Dec. 2001),
  New York University Stern School of Business.
  
  Engle, R.F. (2002),
  \dQuote{Dynamic Conditional Correlation: A Simple Class of
  Multivariate Generalized Autoregressive Conditional
  Heteroskedasticity Models.}
  \emph{Journal of Business and Economic Statistics}
  \bold{20}, 339--350.
}

\seealso{
\code{\link{summary.ccc}}, 
\code{\link{simulateCCC}}, 
\code{\link{plot.ccc}}, 
\code{\link[Rsolnp]{solnp}}
}

\examples{
\dontrun{
ndim <- 3
nobs <- 3000

## Setting parameters in the DCC part
a <- 0.05
b <- 0.8

## Setting a correlation matrix
R <- diag(0, ndim)
R[lower.tri(R)] <- c(0.8, 0.3, 0.1)
R <- R + t(R)
diag(R) <- 1

## setting parameters in the GARCH part
a0 <- c(0.05, 0.07, 0.1)
A <- matrix(c(0.06, 0.0, 0.0, 0.07, 0.08, 0.002, 0.003, 0.0, 0.06), 
    ndim, ndim)
B <- matrix(c(0.75, 0.09, 0.03, 0.001, 0.81, 0.003, 0.001, 0.008, 0.84), 
    ndim, ndim)

## Simulating data
sim.data <- simulateCCC(R, a0, diag(diag(A)), diag(diag(B)), nobs)

## Estimating a CCC model
garch_ccc <- estimateCCC(inia=a0, iniA=diag(diag(A)), iniB=diag(diag(B)), 
    data=sim.data$eps, model="diagonal")

## Summarizing the results
summary_garch_ccc <- summary(garch_ccc)

## Plotting items
plot(garch_ccc, item="correlation")   # this returns an error since correlation is constant
plot(garch_ccc, item="volatility")    # for volatility (square root of conditional variance)
plot(garch_ccc, item="std.residuals") # for satandardized residuals
plot(garch_ccc, item="return")        # for return (original data)
}
}

\keyword{ccgarch2}
\keyword{CCC}
\keyword{nonlinear}

