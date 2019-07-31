\name{Model.BCBinProb}
\alias{Model.BCBinProb}
\title{
Probabilities for beta and correlated binomial distributions given p's and scale-factors. 
}
\description{
Calculates the probabilities for beta and correlated binomials given values for p's and 
scale-factors. 
}
\usage{
Model.BCBinProb(parameter, model.type, model.name, link, ntrials, covariates.matrix.p, 
covariates.matrix.scalef = matrix(c(rep(1, nrow(covariates.matrix.p))), ncol = 1), 
offset.p = c(rep(0, length(ntrials))), offset.scalef = c(rep(0, length(ntrials))))
}
\arguments{
  \item{parameter}{
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{model.type}{
Takes one of two values i.e. 'p only' or 'p and scale-factor'. The 'p only' value fits a 
linear predictor function to the parameter p which is the 'm(1)' in equation (6) of 
Faddy and Smith (2012) divided by 'N'. The 'p and scale-factor' value fits linear predictor 
functions to both p and the scale-factor.
}
  \item{model.name}{
The model being fitted is one of the two 'beta binomial' or 'correlated binomial'. 
}
  \item{link}{
Takes one of nine values i.e., 'logit', 'probit', 'cloglog', 'cauchit', 'log',
'loglog', 'double exponential', 'double reciprocal', 'power logit'. 
The default is 'cloglog'. The 'power logit' has an attribute of 'power' for
which the default is 1 i.e., a logit link.
}
  \item{ntrials}{
This is a scalar representing the denominator i.e., the length of the probability 
mass function returned is this scalar + 1.
}
  \item{covariates.matrix.p}{
A matrix of covariates for p where rows are the number of values in listbinary and 
columns the covariates. This matrix is extracted from the formulae in function BinaryEPPM. 
However, in the accompanying example it is shown how it can be constructed independently 
of function BinaryEPPM. 
}
  \item{covariates.matrix.scalef}{
A matrix of covariates for the scale-factor where rows are the number of values in listbinary 
and columns the covariates. The default is a vector of ones. This matrix is extracted from the 
formulae in function BinaryEPPM. However, in the accompanying example it is shown how it can 
be constructed independently of function BinaryEPPM. 
}
  \item{offset.p}{
An offset vector for p. The default is a vector of ones.
}
  \item{offset.scalef}{
An offset vector for the scale-factor. The default is a vector of ones.
}
}
\value{
List of arguments input together with a list of probabilities vectors and a data frame of 
values of p, theta (beta binomial) or rho (correlated binomial) and the limits for theta or rho. 
  \item{model}{
The model is either 'beta binomial' or 'correlated binomial'. 
}
  \item{link}{
The link is either 'logit' or 'cloglog'.
}
  \item{parameter}{
A vector of the parameters of the model which is set to initial estimates on function call.
}
  \item{probabilities}{
A list of the vectors of probabilities of the model.
}
  \item{probabilities}{
A data frame of values of p, theta (beta binomial) or rho (correlated binomial) and 
the limits for theta or rho. 
}
}
\references{
Hughes G, Madden L (1995). Some methods allowing for aggregated patterns of disease 
incidence in the analysis of data from designed experiments. Plant Pathology, 44, 927-943.

Kupper L, Haseman J (1978). The use of a correlated binomial model for the analysis 
of toxicological epxeriments. Biometrics, 34(1), 69-76.

Williams D (1996). Overdispersion in logistic linear models. In B Morgan (ed.), 
Statistics in Toxicology, 75-84. Oxford Science Publications.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
link <- 'cloglog'
attr(link, which="p") <- make.link(link)
parameter <- c(-0.68294630,0.03451481)
names(parameter) <- c('p','rho')
model.type <- 'p and scale-factor'
model.name <- 'correlated binomial'
ntrials    <- list(c(rep(10,11)))
Model.BCBinProb(parameter, model.type, model.name, link, ntrials,
                covariates.matrix.p = matrix(c(1),nrow=1), 
                covariates.matrix.scalef = matrix(c(1),nrow=1),
                offset.p = c(0), offset.scalef = c(0))
}
\keyword{ models }
