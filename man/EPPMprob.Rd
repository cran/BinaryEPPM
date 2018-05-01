\name{EPPMprob}
\alias{EPPMprob}
\title{
Calculation of vector of probabilities for a extended Poisson process model (EPPM).
}
\description{
Calculates a vector of probabilities given a vector of rates using the matrix 
exponential function from Goulet, Dutang, Maechler, Firth, Shapira, Stadelmann (2014) 
which depends on the package Matrix of Bates and Maechler (2012).
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
Bates D, Maechler M (2016). Matrix: Sparse and Dense Matrix Classes and Methods. 
R package version 1.2-4, \url{https://CRAN.R-project.org/package=Matrix}.  

Goulet V, Dutang C, Maechler M, Firth D, Shapira M, Stadelmann M. (2014).
Package 'expm'. R package version 0.99-1.1, 
\url{https://CRAN.R-project.org/package=expm}. 

Smith D, Faddy M (2014). CountsEPPM: Mean and Variance Modeling of Count Data. 
R package version 2.0, \url{https://CRAN.R-project.org/package=CountsEPPM}. 
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\keyword{ distribution }
