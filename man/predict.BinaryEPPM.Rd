\name{predict.BinaryEPPM}
\alias{predict.BinaryEPPM}
\title{
Prediction Method for BinaryEPPM Objects
}
\description{
Extract various types of predictions from BinaryEPPM regression models.
}
\usage{
\method{predict}{BinaryEPPM}(object, newdata = NULL, type = c("response", 
     "linear.predictor.p", "linear.predictor.scale.factor",
 "p", "scale.factor", "scale.factor.limits", "mean", 
 "variance",   "distribution", "distribution.parameters"), na.action = na.pass, \dots)
}
\arguments{
  \item{object}{
fitted model object of class "BinaryEPPM".
}
  \item{newdata}{
optionally, a data frame in which to look for variables with which to predict. If omitted, 
the original observations are used.
}
  \item{type}{
character indicating type of predictions: fitted means of responses ("response"), 
   linear predictors ("linear.predictor.p", "linear.predictor.scale.factor"), 
   fitted value of probability of success ("p"),
   fitted value of scale-factor ("scale.factor"), fitted value of mean ("mean"),
scale factor limits ("scale.factor.limits"), 
   fitted value of variance ("variance"), fitted probability distribution ("distribution"),
   parameters of fitted distributions ("distribution.parameters")
}
  \item{na.action}{
function determining what should be done with missing values in \emph{newdata}. The 
default is to predict NA.
}
  \item{\dots}{
some methods for this generic function require additional arguments.
}
}
\value{
A vector or list of the predicted values from the fitted model object.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\seealso{
\code{\link[betareg]{predict.betareg}} 
}
\examples{
data("ropespores.case")
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = 'p only', model.name = 'binomial')                 
predict(output.fn, type = "response")
predict(output.fn, type = "linear.predictor.p")
}
\keyword{ methods }
\keyword{ distribution }
