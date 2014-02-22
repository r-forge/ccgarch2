\name{residDiag}
\alias{residDiag}
\alias{print.residDiag}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Residual diagnostics
%%  ~~function to do ... ~~
}
\description{
Carrying out residual diagnostics for the standardized residuals from 
the CCC-, DCC- or CDCC-GARCH models.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
residDiag(object)
\method{print}{residDiag}(x, digits = 4, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{an object of the class \code{"ccc"}, \code{"dcc"} or 
    \code{"cdcc"}}
\item{x}{an object of the class "residDiag".}
\item{digits}{the number of digits to be displayed.}
\item{\dots}{optional arguments. Currently not in use.}

}
\details{
A convenient wrapper to \code{lb.test} and \code{jb.test} carrying out 
the Jarque-Bera and Ljung-Box tests for the standardized residuals. 
}
\value{
An S3 class object "residDiag" is returned. The object is a list with 
the following components.
\item{jb}{the Jarque-Bera test statistics and p-values for the standardized residuals.}
\item{lb}{the Ljung-Box test statistics and p-values for the standardized residuals.}
\item{lb2}{the Ljung-Box test statistics and p-values for the squared standardized residuals.}
}
\references{
Jarque, C.M. and A.K. Bera (1987),
\dQuote{A Test for Normality of Observations and Regression Residuals},
\emph{International Statistical Review}, \bold{55}, 163--172.

Ljung, G.M. and G.E.P. Box (1978):
\dQuote{On a Measure of Lack of Fit in Time-Series Models},
\emph{Biometrika}, \bold{65}, 297--303.

Lomnicki, Z.A. (1961),
\dQuote{Tests for Departure from Normality in the Case of Linear Stochastic Processes},
\emph{Metrika}, \bold{4}, 37--62.

McLeod, A.I., and W.K. Li (1983):
\dQuote{Diagnostic checking ARMA time series models using
  squared-residual autocorrelations}, 
\emph{Journal of Time Series Analysis}, \bold{4}, 269--273.
%% ~put references to the literature/web site here ~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{lb.test}}, \code{\link{jb.test}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
residDiag(object)
}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ccgarch2}
