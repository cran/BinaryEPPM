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
David M. Smith <smithdm1@us.ibm.com>
}
\keyword{ IO }