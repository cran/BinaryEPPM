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
BinaryEPPM(formula, data, subset = NULL, na.action = NULL, 
       weights = NULL, model.type = "p and scale-factor", 
       model.name = "generalized binomial", link = "cloglog", 
       initial = NULL, method = "Nelder-Mead", 
       pseudo.r.squared.type = "square of correlation", control = NULL)
}
\arguments{
  \item{formula}{
Formulae for the probability of a success p and scale-factor. The 
object used is from the package \code{\link[Formula]{Formula}} 
of Zeileis and Croissant (2010) which allows multiple 
parts and multiple responses. "formula" should consist of a 
left hand side (lhs) of single response variable and a right hand 
side (rhs) of one or two sets of variables for the linear predictors 
for the mean and (if two sets) the variance. This is as used for 
the R function "glm" and also, for example, as for the package 
"betareg" (Cribari-Neto and Zeileis, 2010). The function identifies 
from the argument data whether a data frame (as for use of "glm") 
or a list has been input. The list should be exactly the same as 
for a data frame except that the response variable is a list of 
vectors of frequency distributions rather than two vectors of 
paired counts of number responding (r) out of number tested as 
for the data frame. The subordinate functions fit models where 
the response variables are "p.obs", or "scalef.obs" according 
to the model type being fitted. The values for these response 
variables are not input as part of "data", they are calculated 
within the function from a list of grouped binary data input.
If the "model.type" is "p only", "formula" consists of a lhs 
of the response variable and a rhs of the terms of the linear
predictor for the mean model. If the "model.type" is "p and 
scale-factor" there are two sets of terms in the rhs of 
"formula" i.e., "p.obs" and "scalef.obs" together with the 
two sets of terms for the linear predictors of p and scale-factor.
}
  \item{data}{
"data" should be either a data frame (as for use of "glm") 
or a list. The list should be exactly the same as for a 
data frame except that the response variable is a list 
of vectors of frequency distributions rather than a 
vector of single counts as for the data frame. Only one list
is allowed within "data" as it is identified as the dependent
variable. If other lists are in "data", for example for use as
weights, they should be removed from "data" prior to calling 
this function. The extracted list can be called using the
"weights" argument to this function. Within the function 
a working list "listcounts" and data frames with components 
such as "p.obs", "scalef.obs", "covariates", "offset.mean", 
"offset.variance" are set up . The component "covariates" 
is a data frame of vectors of covariates in the model. 
The component "listcounts" is a list of vectors of 
frequency distributions, or the single pairs of r/n in 
grouped form if "data" is a data frame. 
}
  \item{subset}{
Subsetting commands.
}
  \item{na.action}{
Action taken for NAs in data.
}
  \item{weights}{
Vector of list of lists of weights.
}
  \item{model.type}{
Takes one of two values i.e. "p only" or "p and scale-factor". 
The "p only" value fits a linear predictor function to the 
parameter \emph{a} in equation (3) of Faddy and Smith (2012). If 
the model type being fitted is binomial, modeling \emph{a} is the 
same as modeling the mean. For the negative binomial the 
mean is \emph{b exp(a)-1)}, \emph{b} also being as in equation (3) 
of Faddy and Smith (2012). The "p and scale-factor" value 
fits linear predictor functions to both the probability 
of a success p and the scale-factor. 
}
  \item{model.name}{
If model.type is "p only" the model being fitted is one 
of the four "binomial", "Faddy distribution", "beta binomial", "correlated binomial". 
If model.type is "p and scale-factor" the model being fitted is 
either "general" i.e. as equations (4) and (6) of Faddy and Smith (2012)
or one of the two "beta binomial", "correlated binomial". 
}
  \item{link}{
Takes one of nine values i.e., 'logit', 'probit', 'cloglog', 'cauchit', 'log',
'loglog', 'double exponential', 'double reciprocal', 'power logit'. 
The default is 'cloglog'. The 'power logit' has an attribute of 'power' for
which the default is 1 i.e., a logit link.
}
  \item{initial}{
This is a vector of initial values for the parameters. 
If this vector is NULL then initial values based on a 
fitting binomial models using "glm" are calculated 
within the function.
}
  \item{method}{
Takes one of the two values "Nelder-Mead" or "BFGS" these 
being arguments of \code{optim}.
}
  \item{pseudo.r.squared.type}{
Takes one of the three values "square of correlation", "R square" or "max-rescaled R square".
The "default" is as used in Cribari-Neto and Zeileis (2010) and is the square of the 
correlation between the observed and predicted values on the GLM linear predictor scale.
The other two are as described in Cox and Snell (1989), and Nagelkerke (1991) and apply to 
logistic regression.
}
  \item{control}{
"control" is a list of control parameters as used in 
"optim". If this list is NULL the defaults for 
"optim" are set as 
"control <- list(fnscale=-1, trace=0, maxit=1000)". 
The control parameters that can be changed 
by inputting a variable length list are 
"fnscale, trace, maxit, abstol, reltol, alpha, beta, gamma". 
Details of "optim" and its control parameters are 
available in the online R help manuals.
}
}
\value{
\item{data.type}{The type of the data i.e., data frame or list}
\item{list.data}{Data as a list of lists of frequency distributions}
\item{call}{The call of the function}
\item{formula}{The formula argument}
\item{model.type}{The type of model being fitted}
\item{model.name}{The model being fitted}
\item{link}{The link function}
\item{covariates.matrix.p}{The design matrix for the probability of a success}
\item{covariates.matrix.scalef}{The design matrix for the scalefactor}
\item{offset.p}{The offset vector for the probability of a success}
\item{offset.scalef}{The offset vector for the scalefactor}
\item{coefficients}{Estimates of model parameters}
\item{loglikelihood}{Loglikelihood}
\item{vcov}{The variance/covariance matrix}
\item{n}{The number of observations}
\item{nobs}{The number of observations}
\item{df.null}{The degrees of freedom of the null model}
\item{df.residual}{The degrees of freedom of the residual}
\item{vnmax}{Vector of maximums of grouped count data vectors in list.counts}
\item{weights}{Vector or list of weights}
\item{converged}{Whether the iterative process converged, TRUE or FALSE}
\item{iterations}{Number of iterations taken}
\item{method}{Method for optim either Nelder-Mead or BFGS}
\item{pseudo.r.squared}{Pseudo R**2 value}
\item{start}{Starting values for iterative process}
\item{optim}{Estimates of model parameters}
\item{control}{Control parameters for \code{optim}}
\item{fitted.values}{Fitted values for probability of success}
\item{y}{Dependent variable}
\item{terms}{Terms in model fitted}
}
\references{
Cox DR, Snell EJ. (1989). \emph{Analysis of Binary Data}. Second Edition. Chapman \& Hall.

Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.

Grun B, Kosmidis I, Zeileis A. (2012). Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
\emph{Journal of Statistical Software}, \bold{48}(11), 1-25. \doi{10.18637/jss.v048.i11}.

Faddy M, Smith D. (2012). Extended Poisson Process Modeling and 
Analysis of Grouped Binary Data. \emph{Biometrical Journal}, \bold{54}, 426-435.
\doi{10.1002/bimj.201100214}.

Nagelkerke NJD. (1991). A Note on a General Definition of the Coefficient of Determination.
\emph{Biometrika}, \bold{78}, 691-692.

Smith D, Faddy M. (2019). Mean and Variance Modeling of Under-Dispersed and Over-Dispersed
Grouped Binary Data. \emph{Journal of Statistical Software}, \bold{90}(8), 1-20.
\doi{10.18637/jss.v090.i08}.

Zeileis A, Croissant Y. (2010). Extended Model Formulas in R: Multiple Parts and Multiple Responses.
\emph{Journal of Statistical Software}, \bold{34}(1), 1-13. \doi{10.18637/jss.v034.i01}.
}
\seealso{
\code{\link[CountsEPPM:CountsEPPM-package]{CountsEPPM}} 
\code{\link[betareg]{betareg}} 
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\examples{
data("ropespores.case") 
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = "p only", model.name = "binomial")   
summary(output.fn)
}
\keyword{ models }
