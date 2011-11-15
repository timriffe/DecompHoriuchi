DecompContinuous <-
function(func,rates1,rates2,N,...){
	# number of interval jumps   
	L 			<- nrow(rates1) # number of ages
	P 			<- ncol(rates1) # number of factors
	d 			<- (rates2-rates1)
	deltascale 	<- .5:(N-.5)/N
	effectmat 	<- rates1*0
	for (k in 1:N){ # over intervals
		# this part implements the assumption that the other rates (all ages and factors)
		# are also changing linearly in sync
		ratesprop <- rates1 + d * deltascale[k]
		for (i in 1:P){ # across effects
			deltaiak     	<- 0
			deltaia     	<- rep(0,L)
			for (a in 1:L){ # down ages
				# now, we select a single element, [a,i] and increment it forward by .5 delta
				# and also bring it backward by .5 delta
				ratesminus        	<-     ratesplus        <- ratesprop
				ratesplus[a,i]     	<-     ratesplus[a,i] + (d[a, i] / (2 * N))
				ratesminus[a,i] 	<-     ratesminus[a,i] - (d[a, i] / (2 * N))
				# the difference between the funciton evaluated with these rates is the change in y
				# due to this age (La), factor (Pi) and interval (Nk)
				deltaiak         	<-     func(ratesplus,...) - func(ratesminus,...)
				# for this factor and interval, we define a vector of these mini-deltas over age
				deltaia[a]         	<-     deltaiak
			}
			# and when age is done looping we add it to the effects matrix in the appropriate column
			# the the effects matrix sums over all intervals (split according to N- bigger N = finer)
			effectmat[,i] <- effectmat[,i] + deltaia
		}
	}
	return(effectmat)
}

