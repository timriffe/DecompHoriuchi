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

