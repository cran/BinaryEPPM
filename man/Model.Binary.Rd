\name{Model.Binary}
\alias{Model.Binary}
\title{
Function for obtaining output from distributional models.
}
\description{
Produces output of model, parameters and probabilities from the various models.
}
\usage{
Model.Binary(parameter, model.type, model.name, link, ntrials, covariates.matrix.p, 
covariates.matrix.scalef, offset.p, offset.scalef)
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
If model.type is 'p only' the model being fitted is one of the six 'binomial', 'over-dispersed-one', 
'over-dispersed-two', 'EPPM binomial', 'beta binomial' or 'correlated binomial'. If model.type 
is 'p and scale-factor' the model being fitted is one of the three 'EPPM binomial', 
'beta binomial' or 'correlated binomial'.
}
  \item{link}{
Takes one of nine values i.e., 'logit', 'probit', 'cloglog', 'cauchit', 'log',
'loglog', 'double exponential', 'double reciprocal', 'power logit'. 
The default is 'cloglog'. The 'power logit' has an attribute of 'power' for
which the default is 1 i.e., a logit link.
}
  \item{ntrials}{
This is a scalar representing the denominator i.e., the length of the probability mass 
function returned is this scalar + 1.
}
  \item{covariates.matrix.p}{
A matrix of covariates for p where rows are the number of values in listbinary and columns 
the covariates. This matrix is extracted from the formulae in function BinaryEPPM. However, 
in the accompanying example it is shown how it can be constructed independently of function BinaryEPPM. 
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
The output from either Model.BCBinProb, Model.GB, Model.Binary, Model.JMVGB, or Model.ODB.
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
link <- 'cloglog'
attr(link, which="p") <- make.link(link)
parameter <- c(-0.68294630,0.03451481)
names(parameter) <- c('p','rho')
model.type <- 'p and scale-factor'
model.name <- 'correlated binomial'
ntrials    <- list(c(rep(10,11)))
Model.Binary(parameter, model.type, model.name, link, ntrials,
             covariates.matrix.p = matrix(c(1),nrow=1), 
             covariates.matrix.scalef = matrix(c(1),nrow=1),
             offset.p = c(0), offset.scalef = c(0))
}
\keyword{ models }