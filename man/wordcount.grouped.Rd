\name{wordcount.grouped}
\alias{wordcount.grouped}
\docType{data}
\title{
Number of occurences of an article in five-word and ten-word samples from two authors.
}
\description{
The data are the number of occurences of an article in five-word and ten-word samples from 
Macaulay's 'Essay on Milton' and G.K. Chesterton's essay 'About the workers'.
}
\usage{data("wordcount.grouped")}
\format{
  The format is:
List of 4
 $ author      : Factor w/ 2 levels " Macaulay"," Chesterton": 1 1 2 2
 $ fsize       : Factor w/ 2 levels "5","10": 1 2 1 2
 $ vsize       : num [1:4] 5 10 5 10
 $ number.words:List of 4
  ..$ : num [1:6] 45 49 6 0 0 0
  ..$ : num [1:11] 27 44 26 3 0 0 0 0 0 0 ...
  ..$ : num [1:6] 32 35 3 0 0 0
  ..$ : num [1:11] 14 38 16 2 0 0 0 0 0 0 ...
}
\source{
Bailey, B.J.R. (1990). A model for Function Word Counts. Appl. Statist. 39(1), 107-114.
}
\references{
Sellers, K.F., Swift, A.W., Weems, K.S. (2017). A flexible distribution class for count data.
Journal of Statistical Distributions and Applications 41(12), 2616-2626.
}
\examples{
data(wordcount.grouped)
}
\keyword{datasets}
