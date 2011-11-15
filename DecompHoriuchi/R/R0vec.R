R0vec <-
function(x,pfem){
	# where w = the proportion female
	l <- length(x)/2
	sum(x[1:l]*x[(l+1):(2*l)]*pfem)
}

