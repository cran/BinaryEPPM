\name{fitted.BinaryEPPM}
\alias{fitted.BinaryEPPM}
\title{
Extraction of fitted values from BinaryEPPM Objects
}
\description{
This function is generic. Extract the fitted values from models of class "BinaryEPMM".
}
\usage{
\method{fitted}{BinaryEPPM}(object, \dots)
}
\arguments{
  \item{object}{
fitted model object of class "BinaryEPPM".
}
  \item{\dots}{
currently not used.
}
}
\details{
This function is included so that function lrtest from package lmtest can be used.
}
\value{
An vector of class "numeric" of the fitted values from the object of class "BinaryEPMM".
}
\seealso{
\code{\link[stats:stats-package]{fitted}} 
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\keyword{ methods }