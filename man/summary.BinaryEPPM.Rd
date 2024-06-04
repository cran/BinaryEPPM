\name{summary.BinaryEPPM}
\alias{summary.BinaryEPPM}
\title{
Summary of BinaryEPPM Objects
}
\description{
This function is generic. Summary of objects of class "BinaryEPPM".
}
\usage{
\method{summary}{BinaryEPPM}(object, \dots)
}
\arguments{
  \item{object}{
Fitted model object of class "BinaryEPPM".
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\details{
Similar output to that of \code{\link[stats]{summary.glm}} "summary.glm" and \code{\link[betareg]{summary.betareg}}
Cribari-Neto and Zeileis (2010).
}
\value{An object of class "summaryBinaryEPPM" is constructed. This object has the following attributes.
  \item{data.type}{
Indicator of the type of data either 0 "data.frame" or 1 "list".
}
  \item{call}{
The "call" to the function "BinaryEPPM".
}
  \item{formula}{
The model formula in "call".
}
  \item{model.type}{
The model type in "call".
}
  \item{model.name}{
The model name in "call".
}
  \item{link}{
The link function in "call".
}
  \item{offset.p}{
The vector of offsets for the model for p.
}
  \item{offset.scalef}{
The vector of offsets for the model for scale-factor.
}
  \item{coeff.table.p}{
The coefficients of the fitted model for p.
}
  \item{coeff.table.scalef}{
The coefficients of the fitted model for scale-factor.
}
  \item{loglik}{
The log-likelihood of the fitted model.
}
  \item{n}{
The number of observations. Relabelled duplication of "nobs" needed when calling function "lrtest". 
}
  \item{nobs}{
The number of observations.
}
  \item{df.null}{
The degrees of freedom of the null model.
}
  \item{df.residual}{
The degrees of freedom of the residual model.
}
  \item{vnmax}{
Vector of number of "trials" in each observation.
}
  \item{weights}{
Vector of weights for observation.
}
  \item{converged}{
Indicator of convergence.
}
  \item{method}{
Method used to calculate pseudo.r.squared.
}
  \item{pseudo.r.squared}{
The value of the coefficient of determination r squared.
}
  \item{optim}{
Final model fit.
}
  \item{control}{
Control parameters for optimization function "optim".
}
  \item{fitted.values}{
The fitted values.
}
  \item{y}{
The dependent variable in the model.
}
  \item{terms}{
The terms in the model.
}
  \item{npar}{
The number of parameters in the model.
}
}
\seealso{
\code{\link[betareg]{summary.betareg}} 
\code{\link[BinaryEPPM:BinaryEPPM-package]{print.summaryBinaryEPPM}} 
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. 
\doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\keyword{ IO }