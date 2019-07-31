\name{Titanic.survivors.grouped}
\alias{Titanic.survivors.grouped}
\docType{data}
\title{
Titanic survivors data in frequency distribution form.
}
\description{
These data are from the survival log of the Titanic and consist of the number of 
survivors out of the number of passengers broken down into age, sex and class 
categories. The data are in frequency distribution form i.e., a distribution as 
a list of numbers surviving for each age, sex and class category.
}
\usage{data("Titanic.survivors.grouped")}
\format{
  The format is:
List of 4
 $ age           : Factor w/ 2 levels "child","adult": 1 2 1 2 1 2 1 2 1 2 ...
 $ sex           : Factor w/ 2 levels "female","male": 1 1 2 2 1 1 2 2 1 1 ...
 $ pclass        : Factor w/ 3 levels "1st class","2nd class",..: 1 1 1 1 2 2 2 2 3 3 ...
 $ number.survive:List of 12
  ..$ : num [1:2] 0 1
  ..$ : num [1:145] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:6] 0 0 0 0 0 1
  ..$ : num [1:176] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:14] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:94] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:12] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:169] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:32] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:166] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:49] 0 0 0 0 0 0 0 0 0 0 ...
  ..$ : num [1:463] 0 0 0 0 0 0 0 0 0 0 ...
}
\details{
Hilbe (2011) first models these data as a logistic model, then finding 
that they are overdispersed, models them as count data (number of survivors, survive) 
with offset (log of the number of passengers, cases).
}
\source{
Section 9.5, Example 3, pages 263-268, Hilbe, J. (2011). 
}
\references{
Hilbe, J. (2011). Negative Binomial Regression. Cambridge University Press, second edition.
}
\examples{
data(Titanic.survivors.grouped)
}
\keyword{datasets}
