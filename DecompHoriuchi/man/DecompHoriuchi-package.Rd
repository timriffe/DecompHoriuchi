\name{DecompHoriuchi-package}
\alias{DecompHoriuchi-package}
\alias{DecompHoriuchi}
\docType{package}
\title{
2 functions to approximate the continuous decomposition method proposed by Horiuchi et al (2008).}
\description{
This package contains two functions for continuous decomposition. One, a direct R translation of the Matlab code offered by Horiuchi et al (2008) (\code{DecompContinuousOrig()}), and another a slightly friendlier but less flexible implementation (\code{DecompContinuous}).}
\details{
\tabular{ll}{
Package: \tab DecompHoriuchi\cr
Type: \tab Package\cr
Version: \tab 1.1\cr
Date: \tab 2011-07-21\cr
License: \tab General Public\cr
LazyLoad: \tab yes\cr
}
}
\author{
Tim Riffe

Maintainer: Tim Riffe <tim.riffe@gmail.com>
}
\references{
Horiuchi, Wilmoth and Pletcher (2008) A Decomposition Method Based on a Model of Continuous Change. Demography. Vol. 45, (4) pp 785-801
Also see the supplementary materials at \link{http://www.demog.berkeley.edu/~jrw/Papers/decomp.suppl.pdf}, which includes the MatLab program from the authors.
}
\keyword{ package }
\keyword{ stats }
\examples{
library(DecompHoriuchi)
data(rates1)
data(rates2)
################################################################
# DecompContinuous()
Epsilon <- DecompContinuous(R0,rates1,rates2,N=20,pfem=.49)
# how to interpret?:
# the difference in R0 given these rates is:
R0(rates2,.49)-R0(rates1,.49)
# which is also the SUM of Epsilon:
sum(Epsilon)
# and when we look at the individual components of Epsilon, we see their contributions to the change:
# easiest to call barplot()

# 1) flip the matrix on its side (wide)
Epsilon <- t(Epsilon)
# 2) separate negatives and positive:
EPSpos <- Epsilon * .5*(sign(Epsilon)+1) 		
EPSneg <- Epsilon * .5*abs(sign(Epsilon)-1)
# 3) call barplot
barplot(EPSpos,width=rep(1,ncol(Epsilon)),space=0,ylim=range(Epsilon),main="A fake decomposition of R0",col=c("yellow","green"),axisnames=F,xlim=c(0,90),
		ylab="contrib to change in R0",cex.axis=.7)
barplot(EPSneg,width=rep(1,ncol(Epsilon)),add=T,space=0,col=c("yellow","green"),axes=F,axisnames=F)
segments(seq(from=0,to=90,by=10),0,seq(from=0,to=90,by=10),-.027,lty=2,col="grey")
text(seq(from=0,to=90,by=10),-.027,seq(from=0,to=90,by=10),pos=1,xpd=T)
legend("bottomright",fill=c("yellow","green"),legend=c("contrib from change in Lx","contrib from change in Fx"),title="age specific contrib of changes in Fx and Lx",bg="white")
####################################################
# DecompContinuousOrig()
# put rates into vectors,
rates1 <- c(rates1)
rates2 <- c(rates2)
# here we call R0vec instead of R0:
A <- DecompContinuousOrig(func=R0vec,rates1=rates1,rates2=rates2,N=10,pfem=.49)

# reorder A into a matrix (sideways):
A <- t(matrix(A,ncol=2))
# call barplot() twice, once for positive values and again for negative values
Apos <- A * .5*(sign(A)+1) 		
Aneg <- A * .5*abs(sign(A)-1) 	
barplot(Apos,width=rep(1,length(A)/2),space=0,ylim=range(A),main="A fake decomposition of R0",col=c("yellow","green"),axisnames=F,xlim=c(0,90),
		ylab="contrib to change in R0",cex.axis=.7)
barplot(Aneg,width=rep(1,length(A)/2),add=T,space=0,col=c("yellow","green"),axes=F,axisnames=F)
segments(seq(from=0,to=90,by=10),0,seq(from=0,to=90,by=10),-.027,lty=2,col="grey")
text(seq(from=0,to=90,by=10),-.027,seq(from=0,to=90,by=10),pos=1,xpd=T)
legend("bottomright",fill=c("yellow","green"),legend=c("contrib from change in Lx","contrib from change in Fx"),title="age specific contrib of changes in Fx and Lx",bg="white")
# (same plot)
}
