\name{simulateCCC}
\alias{simulateCCC}


\title{Simulating data from a CCC-GARCH process}
\description{This function simulates data from a CCC-GARCH process}
\usage{
simulateCCC(R, a0, A, B, nobs, ncut=1000)
}

\arguments{
    \item{R}{a constant conditional correlation matrix.}
    \item{a0}{a vector of constants in the GARCH part.}
    \item{A}{a matrix of ARCH parameter.}
    \item{B}{a matrix of GARCH parameter.}
    \item{nobs}{the number of observations to be simulated.}
    \item{ncut}{the number of initial entries to be discarded. Default is 1000.}
}

\details{
If the ARCH and GARCH parameter matrices, \code{A} and \code{B}, 
are set as non-diagonal, the corresponding DGP allows for 
interactions in conditional variances (see Nakatani and 
Ter\"asvirta (2009) for instance). 

For estimating a CCC-GARCH model, \code{\link{estimateCCC}} is 
available.
}

\value{

This function returns a list with the following components.

\item{z}{a matrix of the simulated standardized residuals}
\item{h}{a matrix of the simulated conditional variances.}
\item{eps}{a matrix of the simulated time series with CCC-GARCH errors.}

}


\references{

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

  Nakatani, T. and T. Ter{\"{a}}svirta (2009),
  \dQuote{Testing for Volatility Interactions in the Constant 
  Conditional Correlation GARCH Model}, 
  \emph{Econometrics Journal}, \bold{12}, 147--163.
}

\seealso{
\code{\link{estimateCCC}}, 
\code{\link{summary.ccc}}
}

\examples{
## See examples in "estimateCCC".
}


\keyword{ccgarch2}
\keyword{CCC}
\keyword{nonlinear}

