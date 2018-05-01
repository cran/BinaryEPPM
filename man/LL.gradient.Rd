\name{LL.gradient}
\alias{LL.gradient}
\title{
Function used to calculate the first derivatives of the log likelihood with respect to the model parameters.
}
\description{
Function used to calculate the first derivatives of the log likelihood with respect to the model parameters. 
These are numerical derivatives calculated using the numerical derivative functions of Gilbert and Varadhan (2015). 
}
\usage{
LL.gradient(parameter, model.type, model.name, link, ntrials, nsuccess,
            covariates.matrix.p, covariates.matrix.scalef, 
            offset.p, offset.scalef, weights, grad.method)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{parameter}{ 
A vector of the parameters of the model which is set to initial estimates 
on function call.
}
  \item{model.type}{
Takes one of two values i.e. 'p only' or 'p and scale-factor'. The 'p only' 
value fits linear predictor functions to the probability of a success 'p' 
as in Faddy and Smith (2012). The 'p and scale-factor' value fits linear 
predictor functions to both the 'p' and the scale-factor. The default is
'p and scale-factor'.
}
  \item{model.name}{
If model.type is 'p only' the model being fitted is one of the four 'binomial', 
'generalized binomial', 'beta binomial' or 'correlated binomial'. If model.type
is 'p and scale-factor' the model being fitted is one of the three 
'generalized binomial', 'beta binomial' or 'correlated binomial'. Information 
about these models is given in Faddy and Smith (2012). 
The default is 'generalized binomial'.
}
  \item{link}{
Takes one of nine values i.e., 'logit', 'probit', 'cloglog', 'cauchit', 'log',
'loglog', 'double exponential', 'double reciprocal', 'power logit'. 
The default is 'cloglog'. The 'power logit' has an attribute of 'power' for
which the default is 1 i.e., a logit link.
}
  \item{ntrials}{
A vector length 'n+1' representing the number of trials 'n' i.e., a vector 
with all elements equal to 'n'.
}
  \item{nsuccess}{
A vector representing the frequency distribution of the binomial 
distribution for fixed number of trials 'n'. 
}
  \item{covariates.matrix.p}{
A matrix of covariates for the mean where rows are the number of values in 
list.binary and columns the covariates. This matrix is extracted from the 
formulae in function BinaryEPPM. However, in the accompanying example it 
is shown how it can be constructed independently of function BinaryEPPM. 
}
  \item{covariates.matrix.scalef}{
A matrix of covariates for the variance where rows are the number of 
values in list.binary and columns the covariates. The default is a 
vector of ones. This matrix is extracted from the formulae in function 
BinaryEPPM. However, in the accompanying example it is shown how it 
can be constructed independently of function BinaryEPPM. 
}
  \item{offset.p}{
An offset vector for the probability of success p. The default 
is a vector of ones.
}
  \item{offset.scalef}{
An offset vector for the scale-factor. The default is a 
vector of ones.
}
  \item{weights}{
A vector or list of weights for the modeling of probability of success. The default is a 
vector of ones.
}
  \item{grad.method}{
Numerical method used to calculate gradients when the optimization method for optim 
is BFGS either simple or Richardson. This is the grad.method attribute of argument 
method of BinaryEPPM. The default is simple.
}
}
\value{
A vector of numerical first derivatives.
}
\references{
Gilbert P, Varadhan R. (2015). numDeriv: Accurate Numerical Derivatives. R Package version 2014.2-1, 
\url{https://CRAN.R-project.org/package=numDeriv}.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
\dontrun{
    gradient <- grad(LL.Regression.Binary, x = parameter, model.type = model.type,
        model.name = model.name, link = link, ntrials = ntrials, nsuccess = nsuccess,
        covariates.matrix.p = covariates.matrix.p,
        covariates.matrix.scalef = covariates.matrix.scalef, 
        offset.p = offset.p, offset.scalef = offset.scalef, weights = weights,
        grad.method = "Richardson")
    return(gradient)
  }
}
\keyword{ model }

