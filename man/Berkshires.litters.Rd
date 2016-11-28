\name{Berkshires.litters}
\alias{Berkshires.litters}
\docType{data}
\title{
The data are of the number of male piglets born in litters of varying sizes for
the Berkshire breed of pigs. 
}
\description{ 
The data are arranged as a list of binomial frequency distributions where the 
listing is by litter size which is included both as a variate (vsize) and as 
a factor (fsize)}
\usage{data("Berkshires.litters")}
\format{
  The format is:
List of 3
 $ fsize         : Factor w/ 7 levels " size 5"," size 6",..: 1 2 3 4 5 6 7
 $ vsize         : int [1:7] 5 6 7 8 9 10 11
 $ number.success:List of 7
  ..$ : num [1:6] 8 29 72 65 40 3
  ..$ : num [1:7] 5 22 89 129 74 35 4
  ..$ : num [1:8] 1 25 62 131 136 89 26 5
  ..$ : num [1:9] 1 15 79 179 219 149 71 33 4
  ..$ : num [1:10] 2 6 47 117 172 181 117 40 9 2
  ..$ : num [1:11] 2 1 23 65 131 145 120 61 20 3 ...
  ..$ : num [1:12] 0 3 9 22 53 94 72 54 20 4 ...
}
\source{
Brooks, R.J., James, W.H., Gray, E. (1993). Modelling Sub-Binomial Variation in 
the Frequency of Sex Combinations in Litters of Pigs. Biometrics 47, 403-417. 
\doi{10.2307/2532134}
}
\examples{
data("Berkshires.litters")
print(Berkshires.litters)
}
\keyword{datasets}
