\name{rates2}
\alias{rates2}
\docType{data}
\title{
Fake data generated for example.
}
\description{
The first column, Lx, is a discrete survival function for time point 2. The second column are age specific fertility rates.
}
\usage{data(rates2)}
\format{
  The format is:
 num [1:101, 1:2] 0.986 0.982 0.981 0.98 0.979 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "Lx2" "Fx2"
}

\source{
I made up these rates using the package \code{RateSketch}
}

\examples{
library(DecompHoriuchi)
data(rates2)
# nothing fancy
library(DecompHoriuchi)
data(rates1)
data(rates2)
# compare Lx
plot(rates1[,1],type='l',col="blue")
lines(rates2[,1],col="green")
# compare Fx
plot(rates1[,2],type='l',col="blue")
lines(rates2[,2],col="green")
}
\keyword{datasets}
