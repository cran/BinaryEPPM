\name{vcov.BinaryEPPM}
\alias{vcov.BinaryEPPM}
\title{
Variance/Covariance Matrix for Coefficients
}
\description{
Variance/covariance matrix for coefficients of fitted model.
}
\usage{
\method{vcov}{BinaryEPPM}(object, model = c("full", "p", "scale.factor"), \dots)
}
\arguments{
  \item{object}{
fitted model object of class "BinaryEPPM".
}
  \item{model}{
character indicating variance/covariance matrix for all coefficients to be output: all coefficients ("full"),
variance/covariance matrix for coefficients of probability of success ("p"),
variance/covariance matrix for coefficients of scale-factor ("scale.factor")
}
  \item{\dots}{
}
}
\value{
The variance/covariance matrix of the parameters of the fitted model object.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\seealso{
\code{\link[betareg]{betareg}} 
}
\examples{
data("ropespores.case") 
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = 'p only', model.name = 'binomial')   
vcov(output.fn)
}
\keyword{ models }