\name{Hiroshima.grouped}
\alias{Hiroshima.grouped}
\docType{data}
\title{
Data of chromosome aberrations in survivors of Hiroshima grouped into dose ranges 
and represented as frequency distributions.
}
\description{
The data is of chromosome aberrations amongst survivors of the atomic bombs 
exploded over Japan. This data set is grouped into four frequency distributions. 
The data is from Prentice (1986). 
}
\usage{data("Hiroshima.grouped")}
\format{
  The format is:
List of 5
 $ fdose             : Factor w/ 7 levels "0","1-99","100-199",..: 1 2 3 4 5 6 7
 $ vdose             : num [1:7] 0 38 144 244 347 ...
 $ gz                : num [1:7] -0.7829 -0.5704 0.0219 0.5823 1.1573 ...
 $ gzz               : num [1:7] 0.61293 0.32532 0.00048 0.3391 1.33927 ...
 $ number.aberrations:List of 7
  ..$ : num [1:101] 139 66 35 17 3 2 1 0 0 0 ...
  ..$ : num [1:101] 20 23 6 7 3 2 5 2 0 1 ...
  ..$ : num [1:101] 23 12 20 23 6 12 12 12 5 2 ...
  ..$ : num [1:101] 2 2 5 5 3 14 3 2 4 3 ...
  ..$ : num [1:101] 1 1 1 1 3 3 3 3 3 3 ...
  ..$ : num [1:101] 0 0 0 1 1 3 4 1 0 2 ...
  ..$ : num [1:101] 3 1 2 2 2 0 1 0 1 0 ...
}
\details{
The grouped doses of Prentice (1986) have been standardized to gz and gzz with 
mean of 0, standard deviation of 1 form. This was done to match the doses used 
in the related individual survivor data form of the data set Hiroshima.case.Rd.
}
\source{
Prentice RL (1986). "Binary Regression Using an Extended Beta-binomial 
Distribution, with Discussion of Correlation Induced by Covariate Measurement
Errors." Journal of the American Statistical Association, 81, 321-327.
\doi{10.1080/01621459.1986.10478275}
}
\examples{
data("Hiroshima.grouped")
print(Hiroshima.grouped) 
}
\keyword{datasets}
