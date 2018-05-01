\name{print.BinaryEPPM}
\alias{print.BinaryEPPM}
\title{
Printing of BinaryEPPM Objects
}
\description{
Prints objects of class "BinaryEPPM".
}
\usage{
\method{print}{BinaryEPPM}(x, digits = max(3, getOption("digits") - 3), \dots)
}
\arguments{
  \item{x}{
fitted model object of class "BinaryEPPM".
}
  \item{digits}{
digits of printed output.
}
  \item{\dots}{
not currently used.
}
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
print(output.fn)
}
\keyword{ IO }