\name{wordcount.case}
\alias{wordcount.case}
\docType{data}
\title{
Number of occurences of an article in five-word and ten-word samples from two authors.
}
\description{
The data are the number of occurences of an article in five-word and ten-word samples from 
Macaulay's 'Essay on Milton' and G.K. Chesterton's essay 'About the workers'.
}
\usage{data("wordcount.case")}
\format{
  A data frame with 340 observations on the following 5 variables.
  \describe{
    \item{\code{author}}{a factor with levels \code{ Macaulay} \code{ Chesterton}}
    \item{\code{fsize}}{a factor with levels \code{5} \code{10}}
    \item{\code{vsize}}{a numeric vector}
    \item{\code{number.words}}{a numeric vector}
    \item{\code{number.tested}}{a numeric vector}
  }
}
\source{
Bailey, B.J.R. (1990). A model for Function Word Counts. Appl. Statist. 39(1), 107-114.
}
\references{
Sellers, K.F., Swift, A.W., Weems, K.S. (2017). A flexible distribution class for count data.
Journal of Statistical Distributions and Applications 41(12), 2616-2626.
}
\examples{
data(wordcount.case)
}
\keyword{datasets}
