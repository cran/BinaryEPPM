\name{EPPMprob}
\alias{EPPMprob}
\title{
Calculation of vector of probabilities for a extended Poisson process model (EPPM).
}
\description{
Calculates a vector of probabilities given a vector of rates using the matrix 
exponential function from Maechler, Dutang, Goulet, Bates, Firth (2023).
}
\usage{
EPPMprob(vlambda)
}
\arguments{
  \item{vlambda}{
a vector of rates of an extended Poisson process.
}
}
\details{
This is a similar function to that in Smith and Faddy (2014).
}
\value{
The value returned is a vector of probabilities.
}
\references{
Maechler M, Dutang C, Goulet V, Bates D, Firth D. (2023). expm: Matrix Exponential. 
R package version 0.999-8, \url{https://CRAN.R-project.org/package=expm}. 

Smith D, Faddy M (2014). CountsEPPM: Mean and Variance Modeling of Count Data. 
R package version 2.0, \url{https://CRAN.R-project.org/package=CountsEPPM}. 
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\keyword{ distribution }
