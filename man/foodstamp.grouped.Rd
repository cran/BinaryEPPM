\name{foodstamp.grouped}
\alias{foodstamp.grouped}
\docType{data}
\title{
Participation in the federal food stamp program as a list not a data frame.
}
\description{
The dependent variable is a list of frequency distributions of binary variables 
indicating participation in the federal food stamp program. The independent 
variables are two binary ones i.e., tenancy and supplemental income, and a 
continuous one of the log(monthly income+1). 
}
\usage{data("foodstamp.grouped")}
\format{
  The format is:
List of 5
 $ l.participation:List of 150
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  .. [list output truncated]
 $ tenancy        : Factor w/ 2 levels "0","1": 2 2 2 2 1 2 2 2 1 2 ...
 $ suppl.income   : Factor w/ 2 levels "0","1": 1 1 2 1 1 1 1 1 1 1 ...
 $ income         : int [1:150] 271 287 714 521 0 518 458 1266 350 168 ...
 $ l.weights1     :List of 150
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 0.656
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 0.441
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 0.127
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 0.647
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 0 1
  ..$ : num [1:2] 0 0.556
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  ..$ : num [1:2] 1 0
  .. [list output truncated]
}
\source{
Rousseeuw P, Croux C, Todorov V, Ruckstuhl A, Salibian-Barrera M, 
Verbeke T, Koller M, Maechler M (2016).\emph{\bold{robustbase}: Basic Robust Statistics}. 
R package version 0.92-6, \url{CRAN.R-project.org.package=robustbase}.
}
\references{
Kunsch HR, Stefanski LA, Carroll RJ (1989). Conditionally Unbiased Bounded-Influence
Estimation in General Regression Models, with Applications to Generalized Linear Models.
Journal of the American Statistical Association, \bold{84}(406), 460-466.
\doi{10.1080/01621459.1989.10478791}.
}
\examples{
data(foodstamp.grouped)
print(foodstamp.grouped)
}
\keyword{datasets}
