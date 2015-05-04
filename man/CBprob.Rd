\name{CBprob}
\alias{CBprob}
\title{
Calculation of vector of probabilities for the correlated binomial distribution.
}
\description{
Given a vector of parameters and a scalar of the number of trials the function returns a vector of probabilities.
}
\usage{
CBprob(parameter, nt)
}
\arguments{
  \item{parameter}{
A vector of the parameters of the correlated binomial distribution.
}
  \item{nt}{
The number of trials.
}
}
\value{
Vector of probabilities
}
\references{
Kupper L, Haseman J (1978). The Use of a Correlated Binomial Model for the Analysis of Toxicological Experiments. Biometrics, 34(1), 69-76.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
parameter <- c(0.971242852,0.001465007)
names(parameter) <- c('p','rho')
nt <- 37
output <- CBprob(parameter,nt)
print(output)
}
\keyword{ distribution }
