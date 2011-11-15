R0 <-
function(rates,pfem=.49){
	Lx <- rates[,1]
	Fx <- rates[,2]
	sum(Lx*Fx*pfem)
}

