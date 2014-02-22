\name{jb.test}
\alias{jb.test}
\alias{print.jbtest}
\title{ The Lomnicki-Jarque-Bera Test of normality (JB test) }

\description{
This function performs the Lomnicki-Jarque-Bera Test of normality and 
returns test statistics and associated p-values.
}

\usage{
jb.test(obj, ...)

\method{print}{jbtest}(x, digits = 3, \dots)
}

\arguments{
\item{obj}{a vector or matrix of variables to be tested.}
\item{x}{an object of the class "jbtest".}
\item{digits}{the number of digits to be displayed.}
\item{\dots}{optional arguments. Currently not in use.}
}

\details{
\code{obj} is usually a vector or matrix of (the standardized) residuals 
from a GARCH model. 
}


\value{
An S3 class object "jbtest" is returned. It contains the test statistic and 
associated p-values.
}
\references{
Jarque, C.M. and A.K. Bera (1987),
\dQuote{A Test for Normality of Observations and Regression Residuals},
\emph{International Statistical Review}, \bold{55}, 163--172.

Lomnicki, Z.A. (1961),
\dQuote{Tests for Departure from Normality in the Case of Linear Stochastic Processes},
\emph{Metrika}, \bold{4}, 37--62.
}

\seealso{
\code{\link{lb.test}}, \code{\link{residDiag}}
}

\examples{
# for a vector
x <- rnorm(1000)
jb.test(x)

# for a matrix
X <- matrix(rnorm(10000), 5000,2)
jb.test(X)
}

\keyword{ts}
