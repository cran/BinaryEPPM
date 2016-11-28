\name{foodstamp.case}
\alias{foodstamp.case}
\docType{data}
\title{
Participation in the federal food stamp program.
}
\description{
The dependent variable is a binary variable indicating participation in the 
federal food stamp program. The independent variables are two binary ones i.e., 
tenancy and supplemental income, and a continuous one of the log(monthly income+1).
}
\usage{data("foodstamp.case")}
\format{
  A data frame with 150 observations on the following 6 variables.
  \describe{
    \item{\code{participation}}{a numeric vector}
    \item{\code{tenancy}}{a numeric vector}
    \item{\code{suppl.income}}{a numeric vector}
    \item{\code{income}}{a numeric vector}
    \item{\code{n}}{a numeric vector}
    \item{\code{weights1}}{a numeric vector}
  }
}
\source{
Rousseeuw P, Croux C, Todorov V, Ruckstuhl A, Salibian-Barrera M, 
Verbeke T, Koller M, Maechler M (2016). \emph{\bold{robustbase}: Basic Robust Statistics}. 
R package version 0.92-6, \url{CRAN.R-project.org.package=robustbase}.
}
\references{
Kunsch HR, Stefanski LA, Carroll RJ (1989). Conditionally Unbiased Bounded-Influence
Estimation in General Regression Models, with Applications to Generalized Linear Models.
\emph{Journal of the American Statistical Association}, \bold{84}(406), 460-466.
\doi{10.1080/01621459.1989.10478791}.
}
\examples{
data(foodstamp.case)
print(foodstamp.case)
}
\keyword{datasets}
