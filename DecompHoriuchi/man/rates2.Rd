\docType{data}
\name{rates2}
\alias{rates2}
\title{Fake data generated for example.}
\description{
  The first column, Lx, is a discrete survival function for
  time point 1. The second column are age specific
  fertility rates.
}
\examples{
\dontrun{
library(DecompHoriuchi)
data(rates1)
data(rates2)
# nothing fancy
# compare Lx
plot(rates1[,1],type='l',col="blue")
lines(rates2[,1],col="green")
# compare Fx
plot(rates1[,2],type='l',col="blue")
lines(rates2[,2],col="green") }
}
\keyword{datasets}

