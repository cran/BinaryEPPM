\name{BinaryEPPM.distribution}
\alias{BinaryEPPM.distribution}
\title{
Produce summary distributional information for fitted model.
}
\description{
Constructs a data frame of summary distributional information i.e. means, variances and
(optionally) probabilities from the fitted model information output from BinaryEPPM.
}
\usage{
BinaryEPPM.distribution(output.fn, output.probabilities = "no", output.Dparameters = "no")
}
\arguments{
  \item{output.fn}{
This is the object output by BinaryEPPM.
}
  \item{output.probabilities}{
Text input 'no' or 'yes' stating whether or not the vectors of probabilities for the fitted 
model output of output.fn should be included.
}
  \item{output.Dparameters}{
Text input 'no' or 'yes' stating whether or not the vectors of parameters of the underlying 
distributions for the model output of output.fn should be included.
}
}
\details{
The object output by this function is a data frame consisting of means and variances for the 
fitted model of output.fn together (optionally) with related vectors of probabilities and/or 
the vectors of parameters of the underlying distributions. If the model.type is 
'p only' only the values of p, scale-factor, means and variances calculated from the vectors of 
probabilities are output, their variable names being p.obs, scalef.obs, mean.prob and variance.prob
respectively. If the model.type is 'p and scale-factor' then the values of p and scale-factor calculated from
estimates of the parameters and their associated design matrices are also output, in which case 
their variable names are p.par and scalef.par. It is a feature of this package that when the list
form of data inut is used the lengths of the vectors of counts of number of succcesses in 
list.binary (see BinaryEPPM) are considered to be the number of trials, so with this form
of data input there is no need to input the number of trials. 
}
\value{
If the model.type is 'p only' a data frame consisting of the vectors mean.prob, variance.prob
and totalprob is returned. If the model.type is 'p and scale-factor' these three vectors are
returned in the data frame together with two more mean.par and variance.par. If   
output.probabilities is set to 'yes' a list of probability vectors is added to the data frame.
If output.Dparameters is set to 'yes' a list of model parameters is added to the data frame.
}
\author{
David M. Smith <david.m.smith@truvenhealth.com>
}
\examples{
\dontrun{
output.fn <- BinaryEPPM(number.surviving ~ 1+vdose,data=Williams.litters)
output.fn$estses[[1]] <- c('intercept p','slope p dose','intercept scale-factor')
BinaryEPPM.distribution(output.fn,output.probabilities="yes",
                        output.Dparameters = "yes")
}
}
\keyword{ distribution }
