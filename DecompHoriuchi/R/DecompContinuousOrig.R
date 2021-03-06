#' @title DecompContinuousOrig Numeric approximation of Continuous Decomposition
#' 
#' @description This is an exact R implementation of the decomposition code in Matlab offered by the authors in
#' the supplementary material given here: http://www.demog.berkeley.edu/~jrw/Papers/decomp.suppl.pdf. The 
#' difference between \code{DecompContinuous()} and this function is that \code{DecompContinuousOrig} takes 
#' \code{rates1} and \code{rates2} as single vectors, rather than as matrices, and output is also returned as a
#' vector. This difference makes the function more flexible, but may add a step when writing the function to
#' be decomposed. See examples.  
#' 
#' @param func A function specified by the user. This must be able to take the vectors \code{rates1} or 
#' \code{rates2} as its argument, and to return the value of the function, \code{y}, when evaluated for 
#' these rates. It may also have additional arguments, not to be decomposed.
#' @param rates1 vector of covariates to be passed on as arguments to \code{func()}. Covariates
#' can be in any order, as long as \code{func()} knows what to do with them. \code{rates1} is for time 1 
#' (or population 1).
#' @param rates2 is the same as \code{rates1} but for time/population 2.
#' @param N The number of intervals to integrate over. If \code{rates1} are observations from 2005 and 
#' \code{rates2} are observations from 2006 an \code{N} of 20 would imply a delta of 1/20 of a year for 
#' each integration step. Higher \code{N} provides finer results (a lower total residual), but may take 
#' longer to compute. In general, there are decreasing returns to higher \code{N}
#' @param \dots optional parameters to pass on to \code{func()}.
#' 
#' @details The decomposition works by assuming a linear change in all covariates between time 
#' 1 and time 2 (or population 1 and population 2). At each small time step approaching time 2 
#' (the size of which is the inverse of \code{N}) each covariate is moved forward along its linear 
#' trajectory. One at a time, each covariate (of which there are ages*variables of) is switched out 
#' twice, once for its value at 1/(2N) forward and once for its value at 1/(2N) backward in time. The 
#' difference between \code{func()} evaluated with these two rate matrices is the change in \code{y}
#' attributable to that particular covariate and that particular time step. Summing over all N time 
#' steps, we get the contribution to the difference of each covariate, \code{effectmat}. The sum of 
#' \code{effectmat} should come very close to \code{func(rates2)-func(rates1)}. The error decreases 
#' with larger \code{N}, but there is not much point in having an \code{N} higher than 100, and 20 
#' is usually sufficient. This ought to be able to handle a very wide variety of functions. 
#' 
#' @return returns \code{effectmat}, a matrix of the variable effects that is organized in the same way as 
#' \code{rates1} and \code{rates2}. \code{sum(effectmat)} 
#' ought to approximate \code{func(rates2)-func(rates1)}.
#' 
#' @references Horiuchi, Wilmoth and Pletcher (2008) A Decomposition Method Based on a Model of 
#' Continuous Change. Demography. Vol. 45, (4) pp 785-801
#' 
#' @seealso See Also as \code{\link{DecompContinuous}}
#' 
#' @examples 
#' \dontrun{
#' library(DecompHoriuchi)
#' data(rates1)
#' data(rates2)
#' # this version of the function takes the arguments rates1 and rates2 as vectors
#' rates1 <- c(rates1)
#' rates2 <- c(rates2)
#' # look at the function:
#' R0vec
#' # 2 things to point out:
#' # 1) it has an argument pfem, proportion female of births (1/(1+SRB)), that must be specified, but
#' #    that we don't care about decomposing
#' # 2) x is a single vector. The the inside of the function needs to either refer to parts of it by indexing,
#' #    as done here, or else re-assign x to various objects. In this case x[1:l] is Lx and x[(l+1):(2*l)] 
#' # is Fx...
#' A <- DecompContinuousOrig(func=R0vec,rates1=rates1,rates2=rates2,N=10,pfem=.49)
#' # the output, A, is also a single vector. Each element corresponds to the effect of changes in that 
#' # particular covariate toward the overall change in the function.
#' 
#' # this package does not supply default plotting functions, but one strategy might be the following:
#' 
#' # reorder A into a matrix (sideways):
#' A <- t(matrix(A,ncol=2))
#' # call barplot() twice, once for positive values and again for negative values
#' Apos <- A * .5*(sign(A)+1)      # to get a matrix of just the positive values, and zeros in the other cells
#' Aneg <- A * .5*abs(sign(A)-1)   # a matrix for just the negative values
#' 
#' barplot(Apos,width=rep(1,length(A)/2),space=0,ylim=range(A),main="A fake decomposition of R0",
#' col=c("yellow","green"),axisnames=F,xlim=c(0,90), ylab="contrib to change in R0",cex.axis=.7)
#' barplot(Aneg,width=rep(1,length(A)/2),add=T,space=0,col=c("yellow","green"),axes=F,axisnames=F)
#' segments(seq(from=0,to=90,by=10),0,seq(from=0,to=90,by=10),-.027,lty=2,col="grey")
#' text(seq(from=0,to=90,by=10),-.027,seq(from=0,to=90,by=10),pos=1,xpd=T)
#' legend("bottomright",fill=c("yellow","green"),legend=c("contrib from change in Lx",
#' "contrib from change in Fx"),title="age specific contrib of changes in Fx and Lx",bg="white") }
#' 
#' @export 

DecompContinuousOrig <-
function(func,rates1,rates2,N,...){
	# number of interval jumps   
	y1 			<- func(rates1,...)
	y2 			<- func(rates2,...)
	d 			<- rates2-rates1
	n 			<- length(rates1)
	delta 		<- d/N
	x <- rates1+d*matrix(rep(.5:(N-.5)/N,length(rates1)),byrow=TRUE,ncol=N)
	cc <- matrix(0,nrow=n,ncol=N)
	for (j in 1:N){
		for (i in 1:n){
			z 		<- rep(0,n)
			z[i] 	<- delta[i]/2
			cc[i,j] <- func((x[,j]+z),...)-func((x[,j]-z),...)
		}
	}	
	cat("Error of decomposition (in %)\ne =",100*(sum(cc)/(y2-y1)-1),"\n")
	return(rowSums(cc))
}

