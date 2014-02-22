\name{summary.dcc}
\alias{summary.dcc}
\alias{print.summary.dcc}

\title{
Summarizing DCC-GARCH estimation
}

\description{
These functions are S3 methods for class \code{"dcc"} or \code{"summary.dcc"}.
}

\usage{
\method{summary}{dcc}(object, \dots)

\method{print}{summary.dcc}(x, digits = max(3, getOption("digits") - 1), \dots)
}

\arguments{
    \item{object}{an object of class \code{"dcc"}.}
    \item{x}{an object of class \code{"summary.dcc"}.}
    \item{digits}{a number of digits to display.}
    \item{\dots}{optional arguments. Currently not in use.}
}

\details{
  The function \code{summary.dcc()} computes and returns a list of summary
  statistics of the fitted model in \code{object} of the \code{"dcc"} class. 

  Some of the values are printed up to certain decimal places. Actual values of 
  individual components are displayed separately, for instance, by \code{summary(object)$coefficients}. 
  See the \bold{Value} section for a list of components.
  
  For residual diagnostics, \code{\link{residDiag}} carries out the Jarque-Bera test of 
  normality for the standardized residuals and the Ljung-Box test for serial correlations 
  in the standardized and squared standardized residuals. 
}


\value{
  In addition to those available in the \code{object} of the \code{"dcc"} class, 
  the following list components are added.
    \item{nobs}{the number of observations.}
    \item{mu}{the estimate of the conditional mean.}
    \item{garch.par}{the estimates of the parameters in the GARCH part.}
    \item{dcc.par}{the estimates of the parameters in the dcc part.}
    \item{coef}{the table of all estimates with robust s.e. and associated inferencial statistics.}
    \item{convergence}{Integer codes for the convergence status of the 1st and 2nd stage optimization. 
                    See the help of \code{optim} and \code{constrOptim} for detail.}
    \item{counts}{a matrix containing  the number of calls to the objective function and its (numerical) 
                  gradient during the 1st and 2nd stage optimization. The number of calls to the gradient 
                  in the 2nd stage is "NA" because it does not use it. See the help of \code{optim}
                  and \code{constrOptim} for detail.}
    \item{logLik}{a value of the lig-likelihood function at the estimates.}
    \item{AIC}{Akaike's information criterion.}
    \item{BIC}{Bayesian information criterion.}
    \item{CAIC}{Consistent Akaike's information criterion.}
}

\note{
It is time consuming to implement the summary method for \sQuote{\code{dcc}} class 
object. This is because obtaining numerical derivatives, which are used for computing 
the (robust) standard errors of the parameter estimates, is computationally demanding. 
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
}

\seealso{
\code{\link{estimateDCC}}
}

\keyword{ccgarch2}
\keyword{DCC}
\keyword{methods}
