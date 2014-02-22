\name{plot}
\alias{plot.ccc}
\alias{plot.dcc}
\alias{plot.cdcc}

\title{Plotting method for a \code{"ccc"}, \code{"dcc"} or \code{"cdcc"} class object}

\description{
\code{plot} is a generic function for S3 plotting methods for 
a \code{"ccc"}, \code{"dcc"} or \code{"cdcc"} class object. 
}

\usage{
\method{plot}{ccc}(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
  lwd = NULL, \dots)

\method{plot}{dcc}(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
  lwd = NULL, \dots)

\method{plot}{cdcc}(x, item, main = NULL, sub = NULL, xlab = NULL, ylab = NULL, 
  lwd = NULL, \dots)

}

\arguments{
  \item{x}{an object of class \code{"ccc"}, \code{"dcc"} or \code{"cdcc"}, or 
  a corresponding \code{summary} class object.}
  \item{item}{a character string determining the object to be plotted. One of 
     \code{"correlation"}, \code{"volatility"}, \code{"std.residuals"} and \code{"return"}
     is allowed.}
  \item{main}{the main title of the plot.}
  \item{sub}{the sub-title of the plot. If unspecified, no sub-title is displayed.}
  \item{xlab}{the x label of the plot. If missing, \code{"Time"} is used. 
      Setting \code{xlab = ""} displays no x label.}
  \item{ylab}{the y label of the plot. If missing, \code{ylab = "Series 1"} and so on is 
      used. Setting \code{ylab = ""} displays no y label.}
  \item{lwd}{the line width for the plot. If missing, \code{lwd = 1} is used.}
  \item{\dots}{optional arguments passed to \code{plot.zoo()}}
}

\details{
An S3 method \code{plot} is used to plot the return data, estimated 
volatilities or estimated dynamic conditional correlations.

\code{plot} depicts multiple time series plot in a single plot window. 
The argument \code{item} determines which 
item to be plotted. One of \dQuote{\code{correlation}}, 
\dQuote{\code{volatility}}, \dQuote{\code{std.residuals}} and \dQuote{\code{return}} 
should be chosen. 

These are internally coerced to \pkg{zoo} objects, so that optional 
arguments for \code{plot.zoo()} can be passed through the \dQuote{\dots} 
argument.

If \dQuote{\code{correlation}} is selected for a \dQuote{\code{ccc}} class object, 
it is terminated with error.
}

\seealso{
\code{\link{estimateCCC}}, \code{\link{estimateDCC}}, \code{\link{estimateCDCC}}
}

\examples{
## See examples in "estimateCCC", "estimateDCC", or "estimateCDCC".
}


\keyword{ccgarch2}
\keyword{CCC}
\keyword{DCC}
\keyword{CDCC}
