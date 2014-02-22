\name{ccgarch2-package}
\alias{ccgarch2-package}
\docType{package}

\title{ccgarch2: a package for simulation and estimation of conditional correlation 
  GARCH models}
\description{
This package provides functions for simulation and estimation of conditional correlation 
(CC-) GARCH models. The CC-GARCH model includes the Constant Conditional Correlation (CCC-), 
Dynamic Conditional Correlation (DCC-) and corrected Dynamic Conditional Correlation (cDCC-) 
GARCH models.
}

\details{
The package \pkg{ccgarch2} provides functions for simulation and estimation of conditional correlation 
(CC-) GARCH models. The CC-GARCH model includes the Constant Conditional Correlation (CCC-), 
Dynamic Conditional Correlation (DCC-) and corrected Dynamic Conditional Correlation (cDCC-) 
GARCH models.

The CC-GARCH model is defined as follows. 
Let \eqn{\bm{y}_{t}}{y_t} be an \eqn{(N \times 1)}{(N times 1)} vector of returns
\deqn{
  \begin{array}{l}
    \bm{y}_{t} = \bm{\mu} + \bm{\varepsilon}_{t} \\
    \mathsf{E}_{t-1}[\bm{\varepsilon}_{t}] = \bm{0}, \quad
    \mathsf{E}_{t-1}[\bm{\varepsilon}_{t}\bm{\varepsilon}_{t}^{\prime}] = \bm{H}_{t} \\
    \bm{H}_{t} = \bm{D}_{t}^{1/2}\bm{R}_{t}\bm{D}_{t}^{1/2}
  \end{array}
}{ % non-latex environments
  y_t = mu + e_t }
where \eqn{\bm{D}_{t} \equiv \mathsf{diag}(h_{1t}, \ldots, h_{Nt})}{D_{t} = 
diag(h_1t, \ldots, h_Nt)}
is a diagonal matrix whose diagonal entries are the conditional variances of 
\eqn{\bm{y}_{t}}{y_t}, \eqn{\bm{H}_{t}}{H_{t}} is the conditional covariance matrix of 
\eqn{\bm{y}_{t}}{y_t} and \eqn{\bm{R}_{t}}{R_{t}} is the conditional correlation 
matrix of \eqn{\bm{y}_{t}}{y_t} as well as the conditional covariance matrix of 
the standardized residuals, 
\eqn{\bm{z}_{t} \equiv \bm{D}_{t}^{-1/2}\bm{\varepsilon}_{t}}{z_t = D_t^{-1/2}e_t}. 

One of the main properties of the CC-GARCH model 
is to specify the diagonal entries of \eqn{\bm{H}_{t}}{H_{t}} (that is, 
\eqn{h_{it}}{h_it}) and \eqn{\bm{R}_{t}}{R_{t}} separately. While it is 
common to assume that \eqn{h_{it}}{h_it} is described by a variant of 
univariate GARCH models, a number of formulations for \eqn{\bm{R}_{t}}{R_{t}} 
have been proposed. 

In \pkg{ccgarch2}, the GARCH part 
\eqn{\bm{h}_{t}=(h_{1t}, \ldots, h_{Nt})^{\prime}}{h_t = (h_1t, \ldots, h_Nt)'} 
may be modelled as a vector GARCH (1,1) specification such that
\deqn{
  \bm{h}_{t} = \bm{a}_{0} + \bm{A}\bm{\varepsilon}_{t-1}^{(2)} + \bm{B}\bm{h}_{t-1}
}{
h_t = a0 + A*e_{t-1}^{(2)} + B*h_{t-1}
}
where \eqn{\bm{A}}{A} and \eqn{\bm{B}}{B} are \eqn{(N \times N)}{(N times N)} 
matrices, and \eqn{\bm{\varepsilon}_{t-1}^{(2)}}{e_{t-1}^{(2)}} denotes 
the elementwise squared. If \eqn{\bm{A}}{A} and \eqn{\bm{B}}{B} are not 
restricted to diagonal, one may estimate extended types of GARCH models 
that allow for interactions in conditional variances through off-diagonal 
elements of the matrices 
(see Nakatani and Ter\"asvirta (2009), for instance). However, estimation 
of the extended models has difficulty in finding estimates of \eqn{\bm{A}}{A} 
and \eqn{\bm{B}}{B} that satisfy the stationarity conditions of the vector 
GARCH process. Therefore, the codes for estimating CC-GARCH models in 
\pkg{ccgarch2} are restricted to dealing with \eqn{\bm{A}}{A} and 
\eqn{\bm{B}}{B} being diagonal while for simulating data from 
a CC-GARCH process one may set non-zero off-diagonal elements. 

For the correlation part, the following three types of specifications are 
incorporated in \pkg{ccgarch2}. 

\strong{Constant Conditional Correlation; CCC-GARCH} (Bollerslev, 1990)\cr
The simplest way of specifying the conditional correlation part is to 
set the conditional correlation constant over time, namely, 
\deqn{\bm{R}_{t} \equiv \bm{R}}{R_t = R}
For estimation of a CCC-GARCH model \code{\link{estimateCCC}} is available 
while \code{\link{simulateCCC}} simulates the data from a CCC-GARCH process. 

Although the CCC-GARCH model has an advantage in reducing the number of 
parameters to be estimated, the assumption of constant conditional correlation 
is regarded too restrictive. 


\strong{Dynamic Conditional Correlation; DCC-GARCH} (Engle, 2002)\cr
To relax the assumption of constant conditional correlation, Engle (2002) 
proposed the dynamic conditional correlation (DCC-) GARCH model. The 
correlation part of the model is specified as follows. 
\deqn{
  \begin{array}{l}
    \bm{R}_{t} = \bm{Q}_{t}^{*-1/2}\bm{Q}_{t}\bm{Q}_{t}^{*-1/2} \\
    \bm{Q}_{t} = (1-\alpha-\beta)\bm{Q} + \alpha \bm{z}_{t}\bm{z}_{t}^{\prime}
      + \beta \bm{Q}_{t-1}
  \end{array}
}{
 R
}
where \eqn{\bm{Q}_{t}^{*}}{Q_t^{*}} is a diagonal matrix whose entries 
consist of diagonal elements of \eqn{\bm{Q}_{t}}{Q_t}.

For \eqn{\bm{Q}_{t}}{Q_t} to be positive definite, it is assumed that 
\eqn{\bm{Q}}{Q} is positive definite, \eqn{\alpha \geq 0}{alpha >= 0}, 
\eqn{\beta \geq 0}{beta >= 0} and \eqn{\alpha + \beta < 1}{alpha + beta < 1}. 

For estimation of a DCC-GARCH model \code{\link{estimateDCC}} is available 
while \code{\link{simulateDCC}} simulates the data from a DCC-GARCH process. 

Engle (2002) intended to formulate \eqn{\bm{Q}}{Q} as the unconditional second moment of 
\eqn{\bm{z}_{t}}{z_t}. As Aielli (2013) shows, however, this turns 
out not to be the case and interpretation of \eqn{\bm{Q}}{Q} is not immediate. 
Furthermore, Aielli (2013) discusses that 
\eqn{\bm{Q} \neq \mathsf{E}[\bm{z}_{t}\bm{z}_{t}^{\prime}]}{Q != E[z_t z_t]} 
results in inconsistency of the estimators in the DCC-GARCH model in general.

\strong{corrected Dynamic Conditional Correlation; cDCC-GARCH} (Aielli, 2013)\cr
To gain consistency of the estimators in the DCC-GARCH model, Aielli (2013) 
suggested the following modification of \eqn{\bm{Q}_{t}}{Q_t} in the DCC-GARCH model:
\deqn{
  \begin{array}{l}
    \bm{R}_{t} = \bm{Q}_{t}^{*-1/2}\bm{Q}_{t}\bm{Q}_{t}^{*-1/2} \\
    \bm{Q}_{t} = (1-\alpha-\beta)\bm{S} + 
        \alpha \bm{Q}_{t-1}^{*}\bm{z}_{t}\bm{z}_{t}^{\prime}\bm{Q}_{t-1}^{*}
      + \beta \bm{Q}_{t-1}
  \end{array}
}{
 R
}
where \eqn{\bm{S}}{S} is a unit-diagonal matrix. With this formulation, 
\eqn{\bm{S} = \mathsf{E}[\bm{z}_{t}\bm{z}_{t}^{\prime}]}{S = E[z_t z_t]} 
holds, so that the consistency of the estimators is regained (see, 
Aielli (2013) for detail). 

For estimation of a cDCC-GARCH model \code{\link{estimateCDCC}} is available 
while \code{\link{simulateCDCC}} simulates the data from a cDCC-GARCH process. 


The function \code{estimateAAA}, where \code{AAA} is replaced with one of \code{CCC}, 
\code{DCC} and \code{CDCC}, returns an S3 class object \code{"ccc"}, \code{"dcc"} 
or \code{"cdcc"}. The output coisists of the outcomes from the first and second 
stage optimization, original data and matrices of standardized residuals, 
conditional variances and conditional correlations. All of them are combined 
as components of a list object. 

Corresponding S3 methods \code{summary} and \code{plot} are 
available for each class object. \code{summary} summarizes the model 
estimation with coefficient table including estimates and inferencial statistics, 
information criteria such as AIC and diagnostic statistics of the standardized 
residuals. \code{plot} depicts multiple time series plot in a single plot window. 
\code{plot} has an option argument \code{item} that determines which 
item to be plotted. User may choose one of \code{"correlation"}, 
\code{"volatility"}, \code{"std.residuals"} and \code{"level"}. 
}

\note{
It is time consuming (usually takes five to tem minutes) to implement the summary 
method for \sQuote{\code{ccc}}, \sQuote{\code{dcc}}, \sQuote{\code{cdcc}} 
class objects. This is because obtaining numerical derivatives, which are 
used for computing the standard errors of the parameter estimates, is 
computationally demanding. 

% \href{http://bioconductor.org/}{Bioconductor}.

}



\references{

  Aielli, G.P. (2013), 
  \dQuote{Dynamic Conditional Correlation: On Properties and Estimation.}
  \emph{Journal of Business and Economic Statistics}
  \bold{31}, 282--299.

  Bollerslev, T. (1990),
  \dQuote{Modeling the Coherence in Short-Run Nominal Exchange 
  Rates: A Multivariate Generalized {ARCH} Approach,} 
  \emph{Review of Economics and Statistics}, \bold{72}, 498--505.

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



\keyword{ccgarch2}
\keyword{package}

