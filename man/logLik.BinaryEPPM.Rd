\name{logLik.BinaryEPPM}
\alias{logLik.BinaryEPPM}
\title{
Extract Log-Likelihood
}
\description{
This function is generic. It is a method for extracting the log-likelihood for objects of class "BinaryEPPM".
}
\usage{
\method{logLik}{BinaryEPPM}(object, \dots)
}
\arguments{
  \item{object}{
fitted model object of class "BinaryEPPM".
}
  \item{\dots}{
some methods for this generic function require additional arguments
}
}
\details{
logLik is most commonly used for a model fitted by maximum likelihood as is done here.
}
\value{
The log likelihood value for the fitted model object.
}
\seealso{
\code{\link[betareg]{logLik.betareg}} 
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\keyword{ methods }