\name{residuals.BinaryEPPM}
\alias{residuals.BinaryEPPM}
\title{
Residuals for BinaryEPPM Objects
}
\description{
This function is generic. Extract various types of residuals from objects of class "BinaryEPPM".
}
\usage{
\method{residuals}{BinaryEPPM}(object, type = c("spearson", "deviance", "pearson",
        "response", "likelihood", "sdeviance"), \dots)
}
\arguments{
  \item{object}{
Fitted model object of class "BinaryEPPM".
}
  \item{type}{
Type of residuals wanted i.e., standardized Pearson "spearson", deviance "deviance", 
Pearson "pearson",response "response", likelihood "likelihood", standardized deviance "sdeviance".
}
  \item{\dots}{
Some methods for this generic function require additional arguments.
}
}
\details{
Residuals as Cribari-Neto and Zeileis (2010).
}
\value{
An vector of class "numeric" of residuals of a specified type from the object of class "BinaryEPMM".
}
\seealso{
\code{\link[betareg]{residuals.betareg}} 
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\keyword{ models }