\name{KupperHaseman.case}
\alias{KupperHaseman.case}
\docType{data}
\title{
Kupper and Haseman example data
}
\description{
Data of the number of deaths out of number of implants for pregnant female mice
for two groups each of size 10.
}
\usage{data("KupperHaseman.case")}
\format{
  A data frame with 20 observations on the following 3 variables.
  \describe{
    \item{\code{Group}}{a factor with levels \code{Control} \code{Treated}}
    \item{\code{Number.Deaths}}{a numeric vector}
    \item{\code{Number.Implants}}{a numeric vector}
  }
}
\source{
Kupper L, Haseman J. (1978). The Use of a Correlated Binomial Model for the 
Analysis of Toxicological Experiments. Biometrics, 34(1), 69-76.
}
\examples{
data("KupperHaseman.case")
}
\keyword{ datasets }
