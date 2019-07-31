\name{ropespores.case}
\alias{ropespores.case}
\docType{data}
\title{
Dilution series for the presence of rope spores.
}
\description{
Dilution series where at each dilution of a suspension of potato flour a number of 
samples were examined for the presence of rope spores. These data are in data frame form. 
}
\usage{data("ropespores.case")}
\format{
  A data frame with 10 observations on the following 5 variables.
  \describe{
    \item{\code{vdilution}}{a numeric vector}
    \item{\code{fdilution}}{a factor with levels \code{0.25} \code{0.5} \code{1} \code{2} \code{4} \code{8} \code{16} \code{32} \code{64} \code{128}}
    \item{\code{logdilution}}{a numeric vector}
    \item{\code{number.spores}}{a numeric vector}
    \item{\code{number.tested}}{a numeric vector}
  }
}
\source{
Finney, D.J. (1971). Statistical Methods in Biological Assay. Griffin, London, 2nd edition.
}
\examples{
data("ropespores.case")
}
\keyword{datasets}
