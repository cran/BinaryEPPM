\name{GasolineYield}
\alias{GasolineYield}
\docType{data}
\title{
Data on gasoline yields.
}
\description{
Data on gasoline yields originally from Prater (1956). Used as an example in
Cribari-Neto and Zeileis (2010).
}
\usage{data("GasolineYield")}
\format{
  A data frame with 32 observations on the following 8 variables.
  \describe{
    \item{\code{percent}}{a numeric vector}
    \item{\code{vthousand}}{a numeric vector}
    \item{\code{gravity}}{a numeric vector}
    \item{\code{pressure}}{a numeric vector}
    \item{\code{temp10}}{a numeric vector}
    \item{\code{temp}}{a numeric vector}
    \item{\code{batch}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5} \code{6} \code{7} \code{8} \code{9} \code{10}}
  }
}
\source{
Prater NH (1956). Estimate Gasoline Yields from crudes. 
\emph{Petroleum Refiner}, \bold{35}(5), 236-238.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. 
\doi{10.18637/jss.v034.i02}.
}
\examples{
data(GasolineYield)
print(GasolineYield)
}
\keyword{datasets}
