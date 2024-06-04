\name{LL.Regression.Binary}
\alias{LL.Regression.Binary}
\title{
Function called by optim to calculate the log likelihood from the probabilities and hence perform the fitting of 
regression models to the binary data.
}
\description{
Fits specified regression models to the data.
}
\usage{
LL.Regression.Binary(parameter,model.type,model.name,link,ntrials,nsuccess,
                     covariates.matrix.p,covariates.matrix.scalef,
                     offset.p,offset.scalef,weights,grad.method)
}
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
'EPPM extended binomial', 'beta binomial' or 'correlated binomial'. If model.type
is 'p and scale-factor' the model being fitted is one of the three 
'EPPM extended binomial', 'beta binomial' or 'correlated binomial'. Information 
about these models is given in Faddy and Smith (2012). 
The default is 'EPPM extended binomial'.
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
A vector or list of weights for the modeling of probability of success. The default is a vector of ones.
}
  \item{grad.method}{
Numerical method used to calculate gradients either simple or Richardson. The default 
is simple.
}
}
\value{
The log likelihood is returned.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\references{
Faddy M, Smith D. (2012). Extended Poisson Process Modeling and 
Analysis of Grouped Binary Data. \emph{Biometrical Journal}, \bold{54}, 426-435.
\doi{10.1002/bimj.201100214}.
}
\examples{
link <- 'cloglog'
attr(link, which="p") <- make.link(link)
nsuccess <- list(c(rep(0,5),352,479,530,291,101,17))
ntrials  <- list(c(rep(10,11)))
parameter <- c(0.06363398,-0.47085362)
LL.Regression.Binary(parameter, model.type = "p and scale-factor",
            model.name = "EPPM extended binomial", link, ntrials, nsuccess, 
            covariates.matrix.p = matrix(c(1), nrow=1),
            covariates.matrix.scalef = matrix(c(1), nrow=1),
            offset.p = c(0), offset.scalef = c(0),
            weights = list(c(rep(1,11)))) 
}
\keyword{ model }

