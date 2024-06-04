\name{plot.BinaryEPPM}
\alias{plot.BinaryEPPM}
\title{
Diagnostic Plots for BinaryEPPM Objects
}
\description{
This function is generic. Various types of standard diagnostic plots can be produced, involving various 
types of residuals, influence measures etc. It is a minorly modified version of the generic plot function 
of \pkg{betareg} with details of the displays given in Cribari-Neto and Zeileis (2010). The same six displays and arguments list as in Cribari-Neto and Zeileis (2010) are used. The six displays are
"Residuals vs indices of obs", "Cook's distance plot", "Leverage vs predicted values",
"Residuals vs linear predictor", "Normal Q-Q plot of residuals", "Predicted vs observed values".
}
\usage{
\method{plot}{BinaryEPPM}(x, which = 1:4,
    caption = c("Residuals vs indices of obs.", "Cook's distance plot",
    "Leverage vs predicted values", "Residuals vs linear predictor",
    "Normal Q-Q plot of residuals", "Predicted vs observed values"),
    sub.caption = " ", main = "", 
    ask = prod(par("mfcol"), 1) < length(which) && dev.interactive(),
    ..., type = "spearson")
}
\arguments{
  \item{x}{
fitted model object of class "BinaryEPPM".
}
  \item{which}{
numeric. If a subset of plots is required, specify a subset of the numbers 1:6.
}
  \item{caption}{
character. Captions to appear above the plots.
}
  \item{sub.caption}{
character. Common title-above figures if there are multiple.
}
  \item{main}{
character. Title to each plot in addition to the above caption.
}
  \item{ask}{
logical. If true, the user is asked before each plot.
}
  \item{\dots}{
other parameters to be passed through to plotting functions.
}
  \item{type}{
character indicating type of residual to be used, see residuals.BinaryEPPM.
}
}
\details{
The plot method for BinaryEPPM objects produces various plots of diagnostic plots similar to those produced by \pkg{betareg}. See Cribari-Neto and Zeileis (2010) for further details of the displays of \pkg{betareg}.
}
\references{
Cribari-Neto F, Zeileis A. (2010). Beta Regression in R.
\emph{Journal of Statistical Software}, \bold{34}(2), 1-24. \doi{10.18637/jss.v034.i02}.
}
\author{
David M. Smith <dmccsmith@verizon.net>
}
\value{No return value.}
\seealso{
\code{\link[betareg]{plot.betareg}} 
}
\examples{
data("ropespores.case") 
output.fn <- BinaryEPPM(data = ropespores.case,
                  number.spores / number.tested ~ 1 + offset(logdilution),
                  model.type = 'p only', model.name = 'binomial')  
plot.BinaryEPPM(output.fn, which = 1, type= "sdeviance")
}

\keyword{ hplot }

