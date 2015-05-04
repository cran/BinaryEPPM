\name{BinaryEPPM.summary}
\alias{BinaryEPPM.summary}
\title{
Print out summary results.
}
\description{
Print out summary of the results of fitting EPPM models to binary data.
}
\usage{
BinaryEPPM.summary(output.fn)
}
\arguments{
  \item{output.fn}{
The list of objects output from BinaryEPPM. 
}
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
\dontrun{
data("Williams.litters")
output.fn <- BinaryEPPM(number.surviving ~ 1+vdose,data=Williams.litters)
output.fn$estses[[1]] <- c('intercept p','slope p dose','intercept scale-factor')
BinaryEPPM.summary(output.fn) 
}
}
\keyword{ print }
