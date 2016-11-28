\name{Parkes.litters}
\alias{Parkes.litters}
\docType{data}
\title{
The data are of the number of male piglets born in litters of varying sizes for
the Parkes breed of pigs. 
}
\description{ 
The data are arranged as a list of binomial frequency distributions where the 
listing is by litter size which is included both as a variate (vsize) and as 
a factor (fsize)}
\usage{data("Parkes.litters")}
\format{
  The format is:
List of 3
 $ fsize         : Factor w/ 7 levels " size 5"," size 6",..: 1 2 3 4 5 6 7
 $ vsize         : int [1:7] 5 6 7 8 9 10 11
 $ number.success:List of 7
  ..$ : num [1:6] 2 20 41 35 14 4
  ..$ : num [1:7] 3 16 53 78 53 18 0
  ..$ : num [1:8] 0 21 63 117 104 46 21 2
  ..$ : num [1:9] 1 8 37 81 162 77 30 5 1
  ..$ : num [1:10] 0 2 23 72 101 83 46 12 7 0
  ..$ : num [1:11] 0 7 8 19 79 82 48 24 10 0 ...
  ..$ : num [1:12] 0 1 3 15 15 33 13 12 8 1 ...
}
\source{
Brooks, R.J., James, W.H., Gray, E. (1993). Modelling Sub-Binomial Variation 
in the Frequency of Sex Combinations in Litters of Pigs. Biometrics 47, 403-417.
\doi{10.2307/2532134}
}
\examples{
data("Parkes.litters")
print(Parkes.litters)
}
\keyword{datasets}
