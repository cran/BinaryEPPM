\name{LL.Regression.Binary}
\alias{LL.Regression.Binary}
\title{
Function called by optim, or by F.Regression.Binary called by nlm, to calculate 
the log likelihood from the probabilities and hence perform the fitting of 
regression models to the binary data.
}
\description{
Fits specified regression models to the data.
}
\usage{
LL.Regression.Binary(parameter, model.type, model, link, ntrials, nsuccess, 
                    covariates.matrix.p, covariates.matrix.scalef, 
                    offset.p, offset.scalef)
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
  \item{model}{
If model.type is 'p only' the model being fitted is one of the four 'binomial', 
'generalized binomial', 'beta binomial' or 'correlated binomial'. If model.type
is 'p and scale-factor' the model being fitted is one of the three 
'generalized binomial', 'beta binomial' or 'correlated binomial'. Information 
about these models is given in Faddy and Smith (2012). 
The default is 'generalized binomial'.
}
  \item{link}{
Takes one of two value i.e., 'logit' or 'cloglog' these representing the 
logit and complementary log-log link functions respectively. 
The default is 'cloglog'.
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
}
\value{
The log likelihood is returned.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\references{
Faddy M, Smith D (2012). Extended Poisson Process Modeling and 
Analysis of Grouped Data. Biometrical Journal, 54(3), 426-435.
}
\examples{
nsuccess <- list(c(rep(0,5),352,479,530,291,101,17))
ntrials  <- list(c(rep(10,11)))
parameter <- c(0.06363398,-0.47085362)
covariates.matrix.p      <- matrix(c(1),ncol=1) 
covariates.matrix.scalef <- matrix(c(1),ncol=1)
offset.p      <- matrix(c(0),ncol=1)
offset.scalef <- matrix(c(0),ncol=1)
output <- LL.Regression.Binary(parameter,model.type="p and scale-factor",
               model="generalized binomial",link='cloglog',
               nsuccess=nsuccess,ntrials=ntrials,
               covariates.matrix.p,covariates.matrix.scalef,
               offset.p,offset.scalef) 
print(output)
}
\keyword{ model }

