\name{Yorkshires.litters}
\alias{Yorkshires.litters}
\docType{data}
\title{
The data are of the number of male piglets born in litters of varying sizes for 
the Yorkshire breed of pigs. 
}
\description{ 
The data are arranged as a list of binomial frequency distributions where the 
listing is by litter size which is included both as a variate (vsize) and as a factor (fsize)}
\usage{data("Yorkshires.litters")}
\format{
  The format is:
List of 3
 $ fsize         : Factor w/ 9 levels " size 5"," size 6",..: 1 2 3 4 5 6 7 8 9
 $ vsize         : int [1:9] 5 6 7 8 9 10 11 12 13
 $ number.success:List of 9
  ..$ : num [1:6] 3 22 30 37 13 5
  ..$ : num [1:7] 7 18 44 62 27 17 4
  ..$ : num [1:8] 2 14 25 63 69 41 12 5
  ..$ : num [1:9] 2 15 32 70 127 90 45 18 1
  ..$ : num [1:10] 0 8 33 63 106 115 62 30 11 1
  ..$ : num [1:11] 0 3 20 49 79 119 91 59 23 4 ...
  ..$ : num [1:12] 0 0 7 20 60 94 100 47 31 9 ...
  ..$ : num [1:13] 0 1 6 16 29 52 66 43 34 22 ...
  ..$ : num [1:14] 0 2 2 2 14 19 44 45 22 13 ...
}
\source{
Brooks, R.J., James, W.H., Gray, E. (1993). Modelling Sub-Binomial Variation in 
the Frequency of Sex Combinations in Litters of Pigs. Biometrics 47, 403-417.
\doi{10.2307/2532134}
}
\examples{
data("Yorkshires.litters")
print(Yorkshires.litters)
}
\keyword{datasets}
