\name{waldtest.BinaryEPPM}
\alias{waldtest.BinaryEPPM}
\title{
Wald Test of Nested Models for BinaryEPPM Objects
}
\description{
\code{waldtest} is a generic function for comparisons of nested (generalized) linear models via Wald tests. 
}
\usage{
\method{waldtest}{BinaryEPPM}(object, \dots, vcov = NULL,
   test = c("Chisq", "F"))
}
\arguments{
  \item{object}{
an object of class "BinaryEPPM".
}
  \item{\dots}{
further object specifications passed to methods. See below for details.
}
  \item{vcov}{
a function for estimating the covariance matrix of the regression coefficients. If only two models are compared it can also be the covariance matrix of the more general model. 
}
  \item{test}{
character specifying whether to compute the large sample Chi-squared statistic (with asymptotic Chi-squared distribution) or the finite sample F statistic (with approximate F distribution).
}
}
\details{
\code{waldtest} is a generic function for comparisons of nested (generalized)linear models via Wald tests. It does not have the same functionality as the versions of \pkg{betareg} and \pkg{lmtest} with a reduced list of arguments. With these caveats, more details can be obtained from the \bold{Details} pages of those packages.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.

Zeileis A, Hothorn T. (2002). Diagnostic Checking in Regression Relationships. \emph{R News}, \bold{2}(3), 7-10. \url{https://CRAN.R-project.org/doc/Rnews/}.
}
\value{
An object of class \code{"anova"} which contains the residual degrees of freedom, the difference in degrees of freedom, Wald statistic (either \code{"Chisq"} or \code{"F"}) and corresponding p value.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\seealso{
\code{\link[lmtest:waldtest]{waldtest}} 
\code{\link[betareg]{betareg}} 
}
\examples{
\dontrun{
data("ropespores.case") 
output.fn <- BinaryEPPM(data = ropespores.case,
  number.spores / number.tested ~ 1 + offset(logdilution),
  model.type = 'p only', model.name = 'binomial')  
output.fn.one <- BinaryEPPM(data = ropespores.case,
  number.spores / number.tested ~ 1 + offset(logdilution),
  model.type = 'p only', model.name = 'beta binomial')   
waldtest(output.fn, output.fn.one, test = c("Chisq", "F"),
  vcov =  vcov)        
}
}
\keyword{ methods }
