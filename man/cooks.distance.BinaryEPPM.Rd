\name{cooks.distance.BinaryEPPM}
\alias{cooks.distance.BinaryEPPM}
\title{
Cook's distance for BinaryEPPM Objects
}
\description{
Calculates Cook's distances for BinaryEPPM objects.
}
\usage{
\method{cooks.distance}{BinaryEPPM}(model, \dots)
}
\arguments{
  \item{model}{
fitted model object of class "BinaryEPPM".
}
  \item{\dots}{
some mthods for this generic function require additional arguments.
}
}
\details{
Cook's distances as in GLMs.
}
\value{
A vector of Cook's distances.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\seealso{
\code{\link[betareg]{cooks.distance.betareg}} 
}
\examples{
data("ropespores.case") 
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = 'p only', model.name = 'binomial')   
cooks.distance(output.fn)
}
\keyword{ methods }

