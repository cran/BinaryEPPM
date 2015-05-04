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
GBprob(parameter, nt)
}
\arguments{
  \item{parameter}{
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
Faddy M, Smith D (2012). Extended Poisson Process Modeling and Analysis of 
Grouped Binary Data. Biometrical Journal, 54(3), 426-435.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
parameter <- c(0.971242852,0.001465007)
names(parameter) <- c('p','b')
nt <- 37
output <- GBprob(parameter,nt)
print(output)
}
\keyword{ distribution }
