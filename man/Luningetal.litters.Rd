\name{Luningetal.litters}
\alias{Luningetal.litters}
\docType{data}
\title{
Number of trials (implantations) in data of Luning, et al., (1966)
}
\description{
The data are arranged as a list of frequency distributions of the number of female mice with 0,1,2, etc., dead where the number of implants 5 to 10, the number of implants being the length of the frequency distribution. Dose and number of implants are both included as both variates (vdose, vlitter.size) and as factors (fdose, flitter.size).
}
\usage{data("Luningetal.litters")}
\format{
  The format is:
List of 3
 $ vdose      : num [1:3] 0 300 600
 $ fdose      : Factor w/ 3 levels "0","300","600": 1 2 3
 $ number.dead:List of 18
  ..$ : num [1:6] 30 27 9 5 0 0
  ..$ : num [1:6] 27 41 32 17 4 0
  ..$ : num [1:6] 16 32 48 49 15 0
  ..$ : num [1:7] 86 51 14 4 1 0 0
  ..$ : num [1:7] 28 47 59 28 6 1 1
  ..$ : num [1:7] 7 35 45 37 20 9 0
  ..$ : num [1:8] 111 73 31 8 1 0 0 0
  ..$ : num [1:8] 31 61 54 20 19 1 0 0
  ..$ : num [1:8] 5 22 27 36 17 9 3 1
  ..$ : num [1:9] 79 44 23 3 0 1 0 0 0
  ..$ : num [1:9] 12 32 24 22 8 1 0 0 0
  ..$ : num [1:9] 1 4 12 11 8 7 0 2 0
  ..$ : num [1:10] 32 29 8 1 0 0 0 0 0 0
  ..$ : num [1:10] 1 6 9 6 1 1 0 0 0 0
  ..$ : num [1:10] 0 0 2 2 2 0 1 0 0 0
  ..$ : num [1:11] 5 5 2 0 0 0 0 0 0 0 ...
  ..$ : num [1:11] 1 2 1 0 0 0 0 0 0 0 ...
  ..$ : num [1:11] 0 0 0 0 0 0 0 1 0 0 ...
}
\source{
Luning K, Sheridan W, Ytterborn K, Gullberg U (1966). The relationship between the number of implantations and the rate of intra-uterine death in mice. Mutation Research, 3, 444-451.
}
\examples{
data("Luningetal.litters")
print(Luningetal.litters)
}
\keyword{datasets}
