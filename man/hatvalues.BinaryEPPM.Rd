\name{hatvalues.BinaryEPPM}
\alias{hatvalues.BinaryEPPM}
\title{
Extraction of hat matrix values from BinaryEPPM Objects
}
\description{
Extract the values of the hat matrix from models of class "BinaryEPMM".
}
\usage{
\method{hatvalues}{BinaryEPPM}(model, \dots)
}
\arguments{
  \item{model}{
fitted model object of class "BinaryEPPM".
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\value{
The calculated hat values for the fitted model. These are used to calculate Cook's distances.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\seealso{
\code{\link[betareg]{hatvalues.betareg}} 
}
\examples{
data("ropespores.case") 
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = 'p only', model.name = 'binomial')   
hatvalues(output.fn)
}
\keyword{ methods }
