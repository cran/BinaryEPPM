\name{ropespores.grouped}
\alias{ropespores.grouped}
\docType{data}
\title{
Dilution series for the presence of rope spores.
}
\description{
Dilution series where at each dilution of a suspension of potato flour a number of samples were examined for the presence of rope spores. These data are in list form.
}
\usage{data("ropespores.grouped")}
\format{
  The format is:
List of 4
 $ vdilution    : num [1:10] 0.25 0.5 1 2 4 8 16 32 64 128
 $ fdilution    : Factor w/ 10 levels "0.25","0.5","1",..: 1 2 3 4 5 6 7 8 9 10
 $ offset.p     : num [1:10] 1.386 0.693 0 -0.693 -1.386 ...
 $ number.spores:List of 10
  ..$ : num [1:6] 0 0 0 0 0 1
  ..$ : num [1:6] 0 0 0 0 0 1
  ..$ : num [1:6] 0 0 0 0 0 1
  ..$ : num [1:6] 0 0 0 0 0 1
  ..$ : num [1:6] 0 0 0 0 1 0
  ..$ : num [1:6] 0 0 0 1 0 0
  ..$ : num [1:6] 0 0 1 0 0 0
  ..$ : num [1:6] 0 0 1 0 0 0
  ..$ : num [1:6] 1 0 0 0 0 0
  ..$ : num [1:6] 1 0 0 0 0 0
}
\source{
Finney, D.J. (1971). Statistical Methods in Biological Assay. Griffin, London, 2nd edition.
}
\examples{
data("ropespores.grouped")
print(ropespores.grouped)
}
\keyword{datasets}
