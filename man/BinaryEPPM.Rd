\name{BinaryEPPM}
\alias{BinaryEPPM}
\title{
Fitting of EPPM models to binary data.
}
\description{
Fits regression models to under- and over-dispersed binary data 
using extended Poisson process models.
}
\usage{
BinaryEPPM(formula, data, model.type = "p and scale-factor", 
       model = "generalized binomial", link = "cloglog", 
       offset = NULL, initial = NULL, 
       optimization.method = "optim", control = NULL)
}
\arguments{
  \item{formula}{
Formulae for the probability of a success p and scale-factor. The 
package 'Formula' of Zeileis and Croissant (2010) which allows multiple 
parts and multiple responses is used. 'formula' should consist of a 
left hand side (lhs) of single response variable and a right hand 
side (rhs) of one or two sets of variables for the linear predictors 
for the mean and (if two sets) the variance. This is as used for 
the R function 'glm' and also, for example, as for the package 
'betareg' (Cribari-Neto and Zeileis, 2010). The function identifies 
from the argument data whether a data frame (as for use of 'glm') 
or a list (as required in Version 1.0 of this function) has been 
input. The list should be exactly the same as for a data frame 
except that the response variable is a list of vectors of frequency 
distributions rather than two vectors of paired counts of number 
responding (r) out of number tested as for the data frame. The 
subordinate functions fit models where the response variables are 
'p.obs', or 'scalef.obs' according to the model type being fitted. 
The values for these response variables are not input as part of 
'data', they are calculated within the function from a list of 
grouped binary data input. If the 'model.type' is 'p only', 
'formula' consists of a lhs of the response variable and a 
rhs of the terms of the linear predictor for the mean model. 
If the 'model.type' is 'p and scale-factor' there are two 
sets of terms in the rhs of 'formula' i.e., 'p.obs' and 
'scalef.obs' together with the two sets of terms for the 
linear predictors of p and scale-factor. If 
'scale.factor.model'='yes' the second response variable 
used by the subordinate functions would be 'scalef.obs'.     
}
  \item{data}{
'data' should be either a data frame (as for use of 'glm') 
or a list. The list should be exactly the same as for a 
data frame except that the response variable is a list 
of vectors of frequency distributions rather than a 
vector of single counts as for the data frame. Within 
the function a working list 'listcounts' and data frames 
with components such as 'p.obs', 'scalef.obs', 'covariates', 
'offset.mean', 'offset.variance' are set up . The component 
'covariates' is a data frame of vectors of covariates in 
the model. The component 'listcounts' is a list of vectors 
of frequency distributions, or the single pairs of r/n in 
grouped form if 'data' is a data frame. 
}
  \item{model.type}{
Takes one of two values i.e. 'p only' or 'p and scale-factor'. 
The 'p only' value fits a linear predictor function to the 
parameter 'a' in equation (3) of Faddy and Smith (2012). If 
the model type being fitted is binomial modeling 'a' is the 
same as modeling the mean. For the negative binomial the 
mean is 'b'(exp('a')-1), 'b' also being as in equation (3) 
of Faddy and Smith (2012). The 'p and scale-factor' value 
fits linear predictor functions to both the probability 
of a success p and the scale-factor. 
}
  \item{model}{
If model.type is 'p only' the model being fitted is one 
of the four 'binomial', 'Faddy distribution', 'beta binomial', 'correlated binomial'. 
If model.type is 'p and scale-factor' the model being fitted is 
either 'general' i.e. as equations (4) and (6) of Faddy and Smith (2012)
or one of the two 'beta binomial', 'correlated binomial'. 
}
  \item{link}{
This is the generalized linear model function linking the 
probability of success p to the linear predictor. It can be 
either logit ('logit') or complementary log-log ('cloglog'). 
The default is 'cloglog'.
}
  \item{offset}{
'offset' is a one or two component object of vectors of 
offsets for p and scale-factor respectively. }
  \item{initial}{
This is a vector of initial values for the parameters. 
If this vector is NULL then initial values based on a 
fitting binomial models using 'glm' are calculated 
within the function.
}
  \item{optimization.method}{
Takes one of the two values 'optim' or 'nlm' these 
being R optimization functions.
}
  \item{control}{
'control' is a list of control parameters as used in 
'optim' or 'nlm'. If this list is NULL the defaults for 
'optim' are set as 
'control <- list(fnscale=-1, trace=0, maxit=1000)'
and for 'nlm' are set as 
'control <- list(fscale=1, print.level=0, stepmax=1, gradtol=1e-8, steptol=1e-10, iterlim=500)'. 
For 'optim' the control parameters that can be changed 
by inputting a variable length list are 
'fnscale, trace, maxit, abstol, reltol, alpha, beta, gamma'. 
For 'nlm' the parameters are 
'fscale, print.level, stepmax, gradtol,steptol, iterlim'. 
Details of 'optim' and 'nlm' and their control parameters are 
available in the online R help manuals.
}
}
\value{
\item{model.type}{The type of model being fitted}
\item{model}{The model being fitted}
\item{covariates.matrix.mean}{The design matrix for the means}
\item{covariates.matrix.variance}{The design matrix for the variances}
\item{offset.mean}{The offset vector for the means}
\item{offset.variance}{The offset vector for the variances}
\item{estimates}{Estimates of model parameters}
\item{vnmax}{Vector of maximums of grouped count data vectors in list.counts}
\item{loglikelihood}{Loglikelihood}
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
Journal of Statistical Software, 34(2), 1-24. URL http://www.jstatsoft.org/v34/i02/.

Grun B, Kosmidis I, Zeileis A. (2012). Extended Beta Regression in R: Shaken, 
Stirred, Mixed, and Partitioned. Journal of Statistical Software, 
48(11), 1-25. URL http://www.jstatsoft.org/v48/i11/.

Faddy M, Smith D. (2012). Extended Poisson process modelling and 
analysis of grouped binary data. Biometrical Journal, 54, 426-435.

Zeileis A, Croissant Y. (2010). Extended Model Formulas in R: 
Multiple Parts and Multiple Responses. Journal of Statistical Software, 
34(XYZ), 1-12. URL http://www.jstatsoft.org/v34/iXYZ/.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
data("ropespores.case") 
logdilution <- ropespores.case$logdilution
initial  <- c(1.9231205,-0.8964434)
output.fn <- BinaryEPPM(number.spores/number.tested ~ 1,initial=initial,
                  data=ropespores.case,offset=list(logdilution))
output.fn$estses[[1]] <- c('intercept p','intercept scale-factor')
BinaryEPPM.summary(output.fn)
}
\keyword{ models }
