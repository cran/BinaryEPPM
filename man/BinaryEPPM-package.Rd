\name{BinaryEPPM-package}
\alias{BinaryEPPM-package}
\docType{package}
\title{
\packageTitle{BinaryEPPM}
}
\description{
Under- and over-dispersed binary data are modeled using an extended Poisson 
process model (EPPM) appropriate for binary data. A feature of the model is that the 
under-dispersion relative to the binomial distribution only needs to be greater than
zero, but the over-dispersion is  restricted compared to other distributional models  
such as the beta and correlated binomials. Because of this, the examples focus on 
under-dispersed data and how, in combination with the beta or correlated distributions,
flexible models can be fitted to data displaying both under- and over-dispersion. Using
Generalized Linear Model (GLM)  terminology, the functions utilize linear predictors for
the probability of success and scale-factor with various link functions for p, and log 
link for scale-factor, to fit a variety of models relevant to areas such as bioassay.
Details of the EPPM are in Faddy and Smith (2012) and Smith and Faddy (2019). Two 
important changes from version 2.3 are the change to scale-factor rather than 
variance modeling, and the inclusion of a vignette.
}
\details{
\packageIndices{BinaryEPPM}
}
\author{
\packageAuthor{BinaryEPPM}

Maintainer: \packageMaintainer{BinaryEPPM}
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.

Faddy M, Smith D. (2012). Extended Poisson Process Modeling and 
Analysis of Grouped Binary Data. \emph{Biometrical Journal}, \bold{54}, 426-435.
\doi{10.1002/bimj.201100214}.

Grun B, Kosmidis I, Zeileis A. (2012). Extended Beta Regression in R: Shaken, Stirred, Mixed, and Partitioned.
\emph{Journal of Statistical Software}, \bold{48}(11), 1-25. \doi{10.18637/jss.v048.i11}.

Smith D, Faddy M. (2019). Mean and Variance Modeling of Under-Dispersed and Over-Dispersed
Grouped Binary Data. \emph{Journal of Statistical Software}, \bold{90}(8), 1-20.
\doi{10.18637/jss.v090.i08}.

Zeileis A, Croissant Y. (2010). Extended Model Formulas in R: Multiple Parts and Multiple Responses.
\emph{Journal of Statistical Software}, \bold{34}(1), 1-13. \doi{10.18637/jss.v034.i01}.
}
\keyword{ package }
\seealso{
\code{\link[CountsEPPM:CountsEPPM-package]{CountsEPPM}} 
\code{\link[betareg]{betareg}} 
}
\examples{
data("ropespores.case")
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = 'p only', model.name = 'binomial')                 
summary(output.fn) 
}


