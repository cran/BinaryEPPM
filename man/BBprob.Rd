\name{BBprob}
\alias{BBprob}
\title{
Calculation of vector of probabilities for the beta binomial distribution.
}
\description{
Given a vector of parameters and a scalar of the number of trials the function 
returns a vector of probabilities.
}
\usage{
BBprob(twoparameter, nt)
}
\arguments{
  \item{twoparameter}{
A vector of the parameters of the beta binomial distribution.
}
  \item{nt}{
The number of trials.
}
}
\value{
Vector of probabilities
}
\references{
Smith D (1982). Algorithm AS189. Maximum Likelihood Estimation of the 
Parameters of the Beta Binomial Distribution. Applied Statistics, 32, 196-204.

Williams D (1996). "Overdispersion in Logistic Linear Models." In B Mrgan (ed.), 
Statistics in Toxicology, pp75-84. Oxford Science Publications.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
twoparameter <- c(0.96477815,0.7561417)
names(twoparameter) <- c('p','theta')
nt <- 37
output <- BBprob(twoparameter,nt)
print(output)
}
\keyword{ distribution }