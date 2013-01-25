\name{DecompContinuous}
\alias{DecompContinuous}
\title{DecompContinuous Numeric approximation of Continuous Decomposition}
\usage{
  DecompContinuous(func, rates1, rates2, N, ...)
}
\arguments{
  \item{func}{A function specified by the user. This must
  take \code{rates1} or \code{rates2} as its argument, and
  return the value of the function, \code{y}, when
  evaluated for these rates. The function can also contain
  other arguments as parameters not to be decomposed.}

  \item{rates1}{These are the covariate matrices to be
  passed on as arguments to \code{func()}. Columns separate
  variables and rows separate ages.}

  \item{rates2}{just like \code{rates1} but for a different
  population, or a different time point.}

  \item{N}{The number of intervals to integrate over. If
  \code{rates1} are observations from 2005 and
  \code{rates2} are observations from 2006 an \code{N} of
  20 would imply a delta of 1/20 of a year for each
  integration step. Higher \code{N} should be
  asymptotically better, but may take longer to compute. In
  general, there are decreasing returns to higher
  \code{N}.}

  \item{...}{Optional items to pass on to \code{func()} but
  not to be decomposed.}
}
\value{
  returns \code{effectmat}, a matrix of the variable
  effects that is organized in the same way as rates1 and
  rates2. Variables across columns and ages down rows.
  \code{sum(effectmat)} ought to approximate
  \code{func(rates2)-func(rates1)}.
}
\description{
  Approximation of the continuous decomposition strategy
  proposed by Horiuchi et. al. (2008) for changes in
  \code{y} over time or some other dimension assuming
  gradual changes in the variables that determine \code{y}
  between the two observation points. This is basically an
  R translation of the example MatLab program provided by
  the authors in the supplementary material to the article,
  which as of this writing, was available at:
  \url{http://www.demog.berkeley.edu/~jrw/Papers/decomp.suppl.pdf}.
  The only major difference is in the way function
  variables must be specified. The MatLab program asks for
  covariates in a single vector. This R implementation asks
  for covariates in a matrix, with columns specifying
  variables, and rows specifying age (or other) classes.
}
\details{
  The decomposition works by assuming a linear change in
  all covariates between time 1 and time 2 (or population 1
  and population 2). At each small time step approaching
  time 2 (the size of which is the inverse of \code{N})
  each covariate is moved forward along its linear
  trajectory. One at a time, each covariate (of which there
  are ages*variables of) is switched out twice, once for
  its value at 1/2N forward and once for its value at 1/2N
  backward in time. The difference between \code{func()}
  evaluated with these two rate matrices is the change in
  \code{y} attributable to that particular covariate and
  that particular time step. Summing over all N time steps,
  we get the contribution to the difference of each
  covariate, \code{effectmat}. The sum of \code{effectmat}
  should come very close to
  \code{func(rates2)-func(rates1)}. The error decreases
  with larger \code{N}, but there is not much point in
  having an \code{N} higher than 100, and 20 is usually
  sufficient. This ought to be able to handle a very wide
  variety of functions.
}
\examples{
\dontrun{
library(DecompHoriuchi)
data(rates1)
data(rates2)
# DecompContinuous()
Epsilon <- DecompContinuous(R0,rates1,rates2,N=20,pfem=.49)
# how to interpret?:
# the difference in R0 given these rates is:
R0(rates2,.49)-R0(rates1,.49)
' # which is also the SUM of Epsilon:
sum(Epsilon)
# and when we look at the individual components of Epsilon, we see their contributions to the change:
# easiest to call barplot()

# 1) flip the matrix on its side (wide)
Epsilon <- t(Epsilon)
# 2) separate negatives and positive:
EPSpos <- Epsilon * .5*(sign(Epsilon)+1)
EPSneg <- Epsilon * .5*abs(sign(Epsilon)-1)
# 3) call barplot
barplot(EPSpos,width=rep(1,ncol(Epsilon)),space=0,ylim=range(Epsilon),main="A fake decomposition of R0",
col=c("yellow","green"),axisnames=F,xlim=c(0,90), ylab="contrib to change in R0",cex.axis=.7)
barplot(EPSneg,width=rep(1,ncol(Epsilon)),add=T,space=0,col=c("yellow","green"),axes=F,axisnames=F)
segments(seq(from=0,to=90,by=10),0,seq(from=0,to=90,by=10),-.027,lty=2,col="grey")
text(seq(from=0,to=90,by=10),-.027,seq(from=0,to=90,by=10),pos=1,xpd=T)
legend("bottomright",fill=c("yellow","green"),legend=c("contrib from change in Lx",
"contrib from change in #' Fx"),title="age specific contrib of changes in Fx and Lx",bg="white") }
}
\references{
  Horiuchi, Wilmoth and Pletcher (2008) A Decomposition
  Method Based on a Model of Continuous Change. Demography.
  Vol. 45, (4) pp 785-801
}
\seealso{
  See Also as \link{DecompContinuousOrig}, which is the
  exact R implementation of the Matlab code provided by
  Horiuchi et. al. (2008), and a bit more flexible.
}

