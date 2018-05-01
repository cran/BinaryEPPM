\name{CBprob}
\alias{CBprob}
\title{
Calculation of vector of probabilities for the correlated binomial distribution.
}
\description{
Given a vector of parameters and a scalar of the number of trials the function 
returns a vector of probabilities.
}
\usage{
CBprob(twoparameter, nt)
}
\arguments{
  \item{twoparameter}{
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
Kupper L, Haseman J (1978). The Use of a Correlated Binomial Model for the 
Analysis of Toxicological Experiments. Biometrics, 34(1), 69-76.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
twoparameter <- c(0.971242852,0.001465007)
names(twoparameter) <- c('p','rho')
nt <- 37
output <- CBprob(twoparameter,nt)
print(output)
}
\keyword{ distribution }
