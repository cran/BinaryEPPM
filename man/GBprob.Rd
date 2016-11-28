\name{GBprob}
\alias{GBprob}
\title{
Calculation of vector of probabilities for the generalized binomial distribution.
}
\description{
Given a vector of parameters and a scalar of the number of trials the function 
returns a vector of probabilities.
}
\usage{
GBprob(twoparameter, nt)
}
\arguments{
  \item{twoparameter}{
A vector of the parameters of the generalized binomial distribution.
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
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
twoparameter <- c(0.971242852,0.001465007)
names(twoparameter) <- c('p','b')
nt <- 37
output <- GBprob(twoparameter,nt)
print(output)
}
\keyword{ distribution }