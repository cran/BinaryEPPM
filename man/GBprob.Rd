\name{GBprob}
\alias{GBprob}
\title{
Calculation of vector of probabilities for the EPPM binomial distribution.
}
\description{
Given a vector of parameters and a scalar of the number of trials the function 
returns a vector of probabilities. The name GBprob is used to avoid confusion
with EPPMprob which is the function calculating the probabilties given the
constructed vector vector of lambdas.
}
\usage{
GBprob(twoparameter, nt)
}
\arguments{
  \item{twoparameter}{
A vector of the parameters of the EPPM binomial distribution.
}
  \item{nt}{
The number of trials.
}
}
\value{
Vector of probabilities
}
\references{
Faddy M, Smith D. (2012). Extended Poisson Process Modeling and 
Analysis of Grouped Binary Data. \emph{Biometrical Journal}, \bold{54}, 426-435.
\doi{10.1002/bimj.201100214}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\examples{
twoparameter <- c(0.971242852,0.001465007)
names(twoparameter) <- c('p','b')
nt <- 37
GBprob(twoparameter,nt)
}
\keyword{ distribution }