\name{simulateDCC}
\alias{simulateDCC}


\title{Simulating data from a DCC-GARCH process}
\description{This function simulates data from a DCC-GARCH process}
\usage{
simulateDCC(a, b, Q, a0, A, B, nobs, ncut=1000)
}

\arguments{
    \item{a}{a constant in the DCC part.}
    \item{b}{a constant in the DCC part.}
    \item{Q}{a covariance matrix of the standardized residuals.}
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

For estimating a DCC-GARCH model, \code{\link{estimateDCC}} is 
available.
}

\value{

This function returns a list with the following components.

\item{DCC}{a matrix of the simulated DCC}
\item{z}{a matrix of the simulated standardized residuals}
\item{Q}{a matrix of the simulated \eqn{Q_{t}}{Q_t}, which is transformed into DCC.}
\item{h}{a matrix of the simulated conditional variances.}
\item{eps}{a matrix of the simulated time series with DCC-GARCH errors.}

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
\code{\link{estimateDCC}}, 
\code{\link{summary.dcc}}
}

\examples{
## See examples in "estimateDCC".
}


\keyword{ccgarch2}
\keyword{DCC}
\keyword{nonlinear}

