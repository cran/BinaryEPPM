\name{herons.case}
\alias{herons.case}
\docType{data}
\title{
Foraging by green-backed herons.
}
\description{
The data are the numbers of successful attempts out of the numbers of attempts at foraging 
by 20 adult and 20 immature green-backed herons. 
}
\usage{data("herons.case")}
\format{
  A data frame with 40 observations on the following 3 variables.
  \describe{
    \item{\code{Group}}{a factor with levels \code{Adult} \code{Immature}}
    \item{\code{Successes}}{a numeric vector}
    \item{\code{Attempts}}{a numeric vector}
  }
}
\source{
Zhu J, Eickhoff J, Kaiser M (2003). Modelling the Dependence between Number of Trials and Success 
Probability in Beta-Binomial-Poisson Mixture Distributions. Biometrics, 59,955-961.
}
\examples{
data("herons.case")
print(herons.case)
}
\keyword{datasets}
