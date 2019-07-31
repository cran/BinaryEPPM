\name{Hiroshima.case}
\alias{Hiroshima.case}
\docType{data}
\title{
Individual case data of chromosome aberrations in survivors of Hiroshima.
}
\description{
The data is of chromosome aberrations amongst survivors of the atomic bombs
exploded over Japan. This data set is of individual survivors. The data is 
from Morel and Neechal (2012). 
}
\usage{data("Hiroshima.case")}
\format{
  A data frame with 648 observations on the following 6 variables.
  \describe{
    \item{\code{m}}{a numeric vector}
    \item{\code{t}}{a numeric vector}
    \item{\code{t65d_gamma}}{a numeric vector}
    \item{\code{t65d_neutron}}{a numeric vector}
    \item{\code{z}}{a numeric vector}
    \item{\code{zz}}{a numeric vector}
  }
}
\details{
The doses vectors z and zz have been standardized to mean of 0, standard deviation of 1 form.
}
\source{
Morel JG, Neechal NK(2012). Overdispersion Models in SAS. SAS Press.
}
\examples{
data("Hiroshima.case")
}
\keyword{datasets}
