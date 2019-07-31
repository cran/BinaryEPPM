\name{Model.JMVGB}
\alias{Model.JMVGB}
\title{
Probabilities for generalized binomial distributions given p's 
and scale-factors. 
}
\description{
Calculates the probabilities for binomial and generalized 
binomial given values for p's and scale-factors. 
}
\usage{
Model.JMVGB(parameter, model.name, link, ntrials, 
            covariates.matrix.p, covariates.matrix.scalef, 
            offset.p = c(rep(0, length(ntrials))), 
            offset.scalef = c(rep(0, length(ntrials))))
}
\arguments{
  \item{parameter}{
A vector of the parameters of the model which is set to initial 
estimates on function call.
}
  \item{model.name}{
The model being fitted is one of the two 'binomial' or 
'generalized binomial'. 
}
  \item{link}{
Takes one of nine values i.e., 'logit', 'probit', 'cloglog', 'cauchit', 'log',
'loglog', 'double exponential', 'double reciprocal', 'power logit'. 
The default is 'cloglog'. The 'power logit' has an attribute of 'power' for
which the default is 1 i.e., a logit link.
}
  \item{ntrials}{
This is a scalar representing the denominator i.e., the 
length of the probability mass function returned is 
this scalar + 1.
}
  \item{covariates.matrix.p}{
A matrix of covariates for p where rows are the number of 
values in listbinary and columns the covariates. This 
matrix is extracted from the formulae in function BinaryEPPM. 
However, in the accompanying example it is shown how it 
can be constructed independently of function BinaryEPPM. 
}
  \item{covariates.matrix.scalef}{
A matrix of covariates for the scale-factor where rows are 
the number of values in listbinary and columns the covariates. 
The default is a vector of ones. This matrix is extracted 
from the formulae in function BinaryEPPM. However, in the 
accompanying example it is shown how it can be constructed 
independently of function BinaryEPPM. 
}
  \item{offset.p}{
An offset vector for p. The default is a vector of ones.
}
  \item{offset.scalef}{
An offset vector for the scale-factor. The default is 
a vector of ones.
}
}
\value{
List of arguments input together with a list of probabilities 
vectors and a data frame of values of a and b of Equation (5) 
of Faddy and Smith (2012). 
  \item{model}{
The model is either 'binomial' or 'generalized binomial'. 
}
  \item{link}{
The link is either 'logit' or 'cloglog'.
}
  \item{parameter}{
A vector of the parameters of the model which is set to 
initial estimates on function call.
}
  \item{probabilities}{
A list of the vectors of probabilities of the model.
}
  \item{Dparameters}{
A data frame of values of a and b of Equation (5) of 
Faddy and Smith (2012).  
}
}
\references{
Faddy M, Smith D. (2012). Extended Poisson Process Modeling and 
Analysis of Grouped Binary Data. \emph{Biometrical Journal}, \bold{54}, 426-435.
\doi{10.1002/bimj.201100214}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
link <- 'cloglog'
attr(link, which="p") <- make.link(link)
parameter <- c(-0.68294630,0.03451481)
names(parameter) <- c('p','scale-factor')
model.name <- 'generalized binomial'
ntrials <- list(c(rep(10,11)))
Model.JMVGB(parameter, model.name, link, ntrials, 
            covariates.matrix.p = matrix(c(1),nrow=1), 
            covariates.matrix.scalef = matrix(c(1),nrow=1), 
            offset.p = c(0), offset.scalef = c(0))
}
\keyword{ models }

