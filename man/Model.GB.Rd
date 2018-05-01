\name{Model.GB}
\alias{Model.GB}
\title{
Probabilities for binomial and generalized binomial distributions given p's and b. 
}
\description{
Calculates the probabilities for binomial and generalized binomial given values 
for p's and b. 
}
\usage{
Model.GB(parameter, model.name, link, ntrials, covariates.matrix.p, 
         offset.p = c(rep(0, length(ntrials))))
}
\arguments{
  \item{parameter}{
A vector of the parameters of the model which is set to initial estimates 
on function call.
}
  \item{model.name}{
The model being fitted is one of the two 'binomial' or 'generalized binomial'. 
}
  \item{link}{
Takes one of nine values i.e., 'logit', 'probit', 'cloglog', 'cauchit', 'log',
'loglog', 'double exponential', 'double reciprocal', 'power logit'. 
The default is 'cloglog'. The 'power logit' has an attribute of 'power' for
which the default is 1 i.e., a logit link.
}
  \item{ntrials}{
This is a scalar representing the denominator i.e., the length of the 
probability mass function returned is this scalar + 1.
}
  \item{covariates.matrix.p}{
A matrix of covariates for p where rows are the number of values in 
listbinary and columns the covariates. This matrix is extracted from 
the formulae in function BinaryEPPM. However, in the accompanying example 
it is shown how it can be constructed independently of function BinaryEPPM. 
}
  \item{offset.p}{
An offset vector for p. The default is a vector of ones.
}
}
\value{
List of arguments input together with a list of probabilities vectors and 
a data frame of values of a and b of Equation (5) of Faddy and Smith (2012). 
  \item{model}{
The model is either 'binomial' or 'generalized binomial'. 
}
  \item{link}{
The link is either 'logit' or 'cloglog'.
}
  \item{parameter}{
A vector of the parameters of the model which is set to initial estimates 
on function call.
}
  \item{probabilities}{
A list of the vectors of probabilities of the model.
}
  \item{Dparameters}{
A data frame of values of a and b of Equation (5) of Faddy and Smith (2012).  
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
parameter <- c(0.9423342,0.5846321)
names(parameter) <- c('p','b')
model.name <- 'generalized binomial'
ntrials <- list(c(rep(10,11)))
output <- Model.GB(parameter, model.name, link, ntrials, 
                   covariates.matrix.p = matrix(c(1),ncol=1), 
                   offset.p = c(0))
print(output)
}
\keyword{ models }


