\name{powerlogit}
\alias{powerlogit}
\title{
Power Logit Link Function
}
\description{
Computes the power logit link function, including its inverse.
}
\usage{
powerlogit(power = 1)
}
\arguments{
  \item{power}{
power value for the power logit link function.
}
}
\value{
The power logit transformation of theta. All logarithms are natural ones, i.e., 
to base e.
}
\references{
Gaudard MA, Karson MJ, Linder E, Tse Sk (1993). Efficient Designs for Estimation 
in the Power Logistic Quantal Response Model." Statistica Sinica, 3, 233-243.
}
\author{
David M. Smith <smithdm1@us.ibm.com>
}
\keyword{ models }
