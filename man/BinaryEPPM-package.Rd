\name{BinaryEPPM-package}
\alias{BinaryEPPM-package}
\title{
Fitting of EPPM models to binary data.
}
\description{
Fits regression models to under- and over-dispersed binary data 
using extended Poisson process models.
}
\details{
\tabular{ll}{
Package: \tab BinaryEPPM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2015-03-06\cr
License: \tab GPL-2 \cr
}
Using Generalized Linear Model (GLM) terminology, the functions 
utilize linear predictors for the probability of success
and scale-factor with logit or complementary log-log link
functions for p, and log link for scale-factor, to fit the 
regression models. 
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\references{
Faddy M, Smith D (2012). Extended Poisson Process Modeling 
and Analysis of Grouped Binary Data. Biometrical Journal, 
54(3), 426-435.

Zeileis A, Croissant Y (2010). Extended Model Formulas in R: 
Multiple Parts and Multiple Responses. Journal of Statistical 
Software, 34(XYZ), 1-12. URL http://www.jstatsoft.org/v34/iXYZ/.
}
\keyword{ package }
\examples{
data("ropespores.grouped") 
logdilution <- - log(ropespores.grouped$vdilution)
initial  <- c(1.866212,9.490613)
output.fn <- BinaryEPPM(number.spores ~ 1,
                    model.type='p only',model='generalized binomial',
                    data=ropespores.grouped,initial=initial,
                    offset=list(logdilution))
output.fn$estses[[1]] <- c('intercept p','intercept scale-factor')
BinaryEPPM.summary(output.fn)
}