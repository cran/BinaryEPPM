\name{coef.BinaryEPPM}
\alias{coef.BinaryEPPM}
\title{
Extraction of model coefficients for BinaryEPPM Objects
}
\description{
Extract the regression model coefficients from models of class "BinaryEPMM".
}
\usage{
\method{coef}{BinaryEPPM}(object, prtpar = c("full", "p", "scale.factor"), \dots)
}
\arguments{
  \item{object}{
fitted model object of class "BinaryEPPM".
}
  \item{prtpar}{
character indicating coefficients of the fitted model to be output:
all coefficients ("full"),
coefficients of the model for probability of success ("p"),
coefficients of the model for scale-factor ("scale.factor")
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\details{
One of a set of standard extractor functions for fitted model objects of class
"BinaryEPPM.
}
\value{
Vector of coefficients of fitted regression model.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\seealso{
\code{\link[betareg]{betareg}} 
}
\examples{
data("ropespores.case")
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution))   
coef(output.fn, prtpar = "full")
coef(output.fn, prtpar = "p")
coef(output.fn, prtpar = "scale.factor")
}
\keyword{ methods }
