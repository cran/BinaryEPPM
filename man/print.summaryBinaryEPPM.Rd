\name{print.summaryBinaryEPPM}
\alias{print.summaryBinaryEPPM}
\title{
Printing of summaryBinaryEPPM Objects 
}
\description{
Prints the objects of class "summaryBinaryEPPM".
}
\usage{
\method{print}{summaryBinaryEPPM}(x, \dots)
}
\arguments{
  \item{x}{
object output by \code{summary.BinaryEPPM}.
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
print(summary(output.fn))
}
\keyword{ IO }