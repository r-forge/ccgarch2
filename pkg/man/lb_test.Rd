\name{lb.test}
\alias{lb.test}
\alias{print.lbtest}
\title{The Ljung-Box Test for Autocorrelations (LB test)}

\description{
This function performs the Ljung-Box test for autocorrelations, 
and returns test statistics and associated p-values. 
}

\usage{
lb.test(obj, ...)

\method{print}{lbtest}(x, digits = 3, \dots)
}

\arguments{
\item{obj}{a matrix of variables to be tested.}
\item{x}{an object of the class "lbtest".}
\item{digits}{the number of digits to be displayed.}
\item{\dots}{optional arguments. Currently not in use.}
}

\details{
\code{obj} is usually a matrix of (the standardized) residuals 
from a GARCH estimation. When \code{obj} is squared residuals, the test is 
equivalent to the McLeord and Li (1983) test. \code{lb.test()} is a 
wrapper function to \code{Box.test()} with \code{type = "Ljung"}. 
}

\value{
An S3 class object "lbtest" is returned. 
It contains the test statistics and associated p-values from lags 5 to 50 by 5.
}

\references{
Ljung, G.M. and G.E.P. Box (1978):
\dQuote{On a Measure of Lack of Fit in Time-Series Models},
\emph{Biometrika}, \bold{65}, 297--303.

McLeod, A.I., and W.K. Li (1983):
\dQuote{Diagnostic checking ARMA time series models using
  squared-residual autocorrelations}, 
\emph{Journal of Time Series Analysis}, \bold{4}, 269--273.
}

\seealso{
\code{\link{jb.test}}, \code{\link{residDiag}}
}

\examples{
X <- matrix(rnorm(10000), 5000,2)
lb.test(X)
}

\keyword{ts}
