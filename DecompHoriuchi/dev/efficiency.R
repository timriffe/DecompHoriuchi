# TODO: Add comment
# 
# Author: triffe
###############################################################################

# where N is how many time intervals
DecompContinuousOrig <-	function(func, rates1, rates2, N, ...){
	# number of interval jumps   
	y1 			<- func(rates1, ...)
	y2 			<- func(rates2, ...)
	d 			<- rates2 - rates1
	n 			<- length(rates1)
	delta 		<- d / N
	x           <- rates1 + d * matrix(rep(.5:(N - .5) / N, length(rates1)), byrow = TRUE, ncol = N)
	cc          <- matrix(0, nrow = n, ncol = N)
	
	for (j in 1:N){
		for (i in 1:n){
			z 		 <- rep(0, n)
			z[i] 	 <- delta[i] / 2
			cc[i, j] <- func((x[, j] + z), ...) - func((x[, j] - z), ...)
		}
	}	
	#cat("Error of decomposition (in %)\ne =", 100*(sum(cc) / (y2 - y1) - 1), "\n")
	return(rowSums(cc))
}

ROvec <- function(x,pfem){
	# where w = the proportion female
	l <- length(x)/2
	sum(x[1:l]*x[(l+1):(2*l)]*pfem)
}

load("/home/triffe/git/DecompHoriuchi/DecompHoriuchi/data/rates1.rda")
load("/home/triffe/git/DecompHoriuchi/DecompHoriuchi/data/rates2.rda")

library(compiler)
DO_c <- cmpfun(DecompContinuousOrig)
ROvec_c <- cmpfun(ROvec)
R0vec(c(rates1),pfem=.48)

install.packages("rbenchmark")

library(rbenchmark)

benchmark(DecompContinuousOrig(R0vec,c(rates1),c(rates2),N=20,pfem=.48),
		DO_c(R0vec,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48),
		DecompContinuousOrig(func=ROvec_c,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48),
		DO_c(func=ROvec_c,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48),
		replications=100
		)
system.time(cc <- DecompContinuousOrig(R0vec,c(rates1),c(rates2),N=20,pfem=.48))
system.time(cc <- DO_c(R0vec,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48))
system.time(cc <- DecompContinuousOrig(func=ROvec_c,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48))
system.time(cc <- DO_c(func=ROvec_c,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48))
system.time(cc <- DCO2(func=ROvec_c,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48))
system.time(cc <- DCO3(func=ROvec_c,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48)

system.time(cc <- DCO4(func = R0vec, rates1 = c(rates1), rates2 = c(rates2), N = 100, pfem = .48, parallel = TRUE))
system.time(cc <- DCO4(func = R0vec, rates1 = c(rates1), rates2 = c(rates2), N = 100, pfem = .48, parallel = FALSE))


system.time(cc <- DCO4.2(func=R0vec,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48))
system.time(cc <- DCO4_c(func=R0vec,rates1=c(rates1),rates2=c(rates2),N=20,pfem=.48))
args(DCO4)

# make 'parallel' a dependency, so it gets installed with package if necessary.
# then create z as a matrix to apply(mar=2) over. j can be a parSapply().
# cc would be the simplified output of parSapply()

# compare this with simply compiling the loops
# also, of course simply compile 'func'hab

DCO2 <-	function(func, rates1, rates2, N, ...){
	# number of interval jumps   
	d 			<- rates2 - rates1
	n 			<- length(rates1)
	delta 		<- d / N
	x           <- rates1 + d * matrix(rep(.5:(N - .5) / N, length(rates1)), byrow = TRUE, ncol = N)
	cc          <- matrix(0, nrow = n, ncol = N)
	
	for (j in 1:N){
		z <- matrix(rep(x[, j], n), ncol = n)
		cc[, j] <- apply(z + diag(delta/2), 2, func, ...) - apply(z - diag(delta/2), 2, func, ...)
	}	
	return(rowSums(cc))
}

DCO3 <-	function(func, rates1, rates2, N, ...){
	# number of interval jumps   
	y1 			<- func(rates1, ...)
	y2 			<- func(rates2, ...)
	d 			<- rates2 - rates1
	n 			<- length(rates1)
	delta 		<- d / N
	x           <- rates1 + d * matrix(rep(.5:(N - .5) / N, length(rates1)), byrow = TRUE, ncol = N)
	#cc          <- matrix(0, nrow = n, ncol = N)
	cc <- sapply(1:N,function(j, x, func, delta, ...){
				z <- matrix(rep(x[, j], n), ncol = n)
		apply(z + diag(delta/2), 2, func, ...) - apply(z - diag(delta/2), 2, func, ...)
			}, x = x, func = func, delta = delta, ...)
#	for (j in 1:N){
#		z <- matrix(rep(x[, j], n), ncol = n)
#		cc[, j] <- apply(z + diag(delta/2), 2, func, ...) - apply(z - diag(delta/2), 2, func, ...)
#	}	
	return(rowSums(cc))
}

DCO4 <-	function(func, rates1, rates2, N, compile.function = TRUE, parallel = TRUE, ...){
	# by default, your specified function will be compiled.
	# 'compiler' NOT installed on the fly.
	if (compile.function && "compiler" %in% installed.packages()){
		func		<- compiler:::cmpfun(match.fun(func))
	}
	# d = change 
	d 			<- rates2 - rates1
	n 			<- length(rates1)
	# assuming linear change:
	delta 		<- d / N
	# and further assuming that all rate components changed linearly, each col of x is a 1/N time step
	rates.d           <- rates1 + d * matrix(rep(.5:(N - .5) / N, length(rates1)), byrow = TRUE, ncol = N)
	d2 			<- diag(delta / 2)
	
	if (parallel && "parallel" %in% installed.packages()){
		cl <- parallel:::makeCluster(parallel:::detectCores())
		cc <- parallel:::parSapply(cl, 1:N, function(j, rates.d, func, d2, ...){
					z <- matrix(rep(rates.d[, j], n), ncol = n)
					apply(z + d2, 2, func, ...) - apply(z - d2, 2, func, ...)
				}, 
				rates.d = rates.d, func = func, d2 = d, ...)
		parallel:::stopCluster(cl)
	} else {
		cc <- sapply(1:N,function(j, rates.d, func, d2, ...){
					z <- matrix(rep(rates.d[, j], n), ncol = n)
					apply(z + d2, 2, func, ...) - apply(z - d2, 2, func, ...)
				}, 
				rates.d = rates.d, func = func, d2 = d2, ...)
	}	
	return(rowSums(cc))
}

DCO4.2 <-	function(func, rates1, rates2, N, ...){
	# by default, your specified function will be compiled.
	# 'compiler' NOT installed on the fly.
	func		<- compiler:::cmpfun(match.fun(func))
	# d = change 
	d 			<- rates2 - rates1
	n 			<- length(rates1)
	# assuming linear change:
	delta 		<- d / N
	# and further assuming that all rate components changed linearly, each col of x is a 1/N time step
	rates.d     <- rates1 + d * matrix(rep(.5:(N - .5) / N, length(rates1)), byrow = TRUE, ncol = N)
	d2 			<- diag(delta / 2)

	require(parallel)
	cl <- makeCluster(parallel:::detectCores())
	cc <- parSapply(cl, 1:N, function(j, rates.d, func, d2, ...){
				z <- matrix(rep(rates.d[, j], n), ncol = n)
				apply(z + d2, 2, func, ...) - apply(z - d2, 2, func, ...)
			}, 
			rates.d = rates.d, func = func, d2 = d, ...)
	stopCluster(cl)

	return(rowSums(cc))
}
DCO4_c <- cmpfun(DCO4)



# silly test: elapsed time over N by parallel status:
#NN <- seq(10,1000,by=10)
#parstat <- matrix(ncol=2,nrow=length(NN))
#rownames(parstat) <- NN
#colnames(parstat) <- c("parallel TRUE", "parallel FALSE")
#
#for (i in 1:length(NN)){
#parstat[i,1] <- mean(replicate(10,system.time(cc <- DCO4(func = R0vec, rates1 = c(rates1), rates2 = c(rates2), N = NN[i], pfem = .48, parallel = TRUE))[3]))
#parstat[i,2] <- mean(replicate(10,system.time(cc <- DCO4(func = R0vec, rates1 = c(rates1), rates2 = c(rates2), N = NN[i], pfem = .48, parallel = FALSE))[3]))
#}
# took like 2 or 3 hrs on my laptop with 2 cores
# results are here:
parstat <- dput(parstat)
structure(c(0.649699999999939, 0.760899999999674, 0.888799999999901, 
				1.04959999999992, 1.21560000000009, 1.32330000000002, 1.4371000000001, 
				1.57740000000031, 1.59400000000023, 1.71619999999966, 1.87759999999998, 
				1.98380000000034, 2.3742000000002, 2.22509999999966, 2.40089999999946, 
				2.59869999999937, 2.72539999999972, 2.86779999999999, 2.88160000000025, 
				2.94720000000052, 3.06870000000017, 3.07410000000018, 3.08310000000019, 
				3.26819999999971, 3.30479999999952, 3.53259999999973, 3.61460000000079, 
				3.78010000000068, 3.9020999999997, 4.01300000000083, 4.12900000000045, 
				4.1752999999997, 4.26889999999949, 4.36419999999962, 4.34350000000049, 
				4.42410000000018, 4.52489999999998, 4.4963000000007, 4.55909999999967, 
				4.62149999999965, 4.73640000000014, 4.93910000000033, 5.1158999999996, 
				5.26290000000081, 5.44359999999979, 5.48160000000062, 5.5706000000002, 
				5.63430000000044, 5.66599999999962, 5.75870000000032, 5.80159999999996, 
				5.9497999999996, 5.94089999999924, 6.02119999999995, 6.16510000000089, 
				6.2781999999992, 6.3143999999993, 6.3656999999992, 6.45120000000024, 
				6.45740000000114, 6.61249999999964, 6.71119999999974, 6.75190000000111, 
				6.74889999999941, 6.9280999999999, 7.0226999999999, 6.99800000000105, 
				7.20680000000029, 7.18079999999936, 7.23129999999946, 7.39709999999977, 
				7.40959999999977, 7.60849999999955, 7.8148999999994, 7.90709999999999, 
				7.90669999999991, 8.02010000000046, 8.15310000000099, 8.11720000000023, 
				8.24629999999961, 8.27719999999972, 8.50210000000043, 8.52389999999978, 
				8.53699999999917, 8.70690000000031, 8.78589999999931, 8.69409999999989, 
				8.84769999999953, 8.93619999999901, 8.96210000000065, 9.16460000000115, 
				9.12839999999924, 9.25839999999989, 9.4929000000011, 9.47870000000039, 
				9.66319999999942, 9.71029999999991, 9.71590000000033, 9.7908000000014, 
				9.84439999999995, 0.143800000000192, 0.274699999999939, 0.415299999999661, 
				0.55049999999992, 0.692799999999625, 0.831599999999889, 0.965500000000247, 
				1.09630000000034, 1.24249999999993, 1.37820000000029, 1.51979999999894, 
				1.78009999999958, 1.78830000000053, 1.93470000000052, 2.05560000000005, 
				2.19830000000002, 2.32609999999913, 2.46470000000008, 2.59739999999947, 
				2.74710000000014, 2.86029999999992, 3.0012999999999, 3.14000000000087, 
				3.29839999999967, 3.40800000000054, 3.5512999999999, 3.69030000000021, 
				3.82250000000022, 3.95139999999956, 4.08290000000015, 4.22250000000058, 
				4.36040000000066, 4.49649999999892, 4.63179999999993, 4.75980000000018, 
				4.89489999999969, 5.03220000000001, 5.18129999999983, 5.31490000000122, 
				5.44840000000077, 5.58559999999925, 5.72169999999933, 5.86079999999965, 
				6.00890000000036, 6.14160000000011, 6.26879999999983, 6.40929999999935, 
				6.55810000000019, 6.68069999999971, 6.8224000000002, 6.97790000000095, 
				7.08899999999958, 7.2380000000001, 7.36820000000007, 7.51390000000029, 
				7.65290000000023, 7.78269999999975, 7.94089999999924, 8.07760000000053, 
				8.19640000000036, 8.34280000000072, 8.48140000000021, 8.61609999999964, 
				8.75249999999978, 8.90139999999992, 9.03429999999935, 9.18679999999949, 
				9.31849999999977, 9.44539999999979, 9.60190000000002, 9.73979999999974, 
				9.8747000000003, 10.0075000000001, 10.1692999999999, 10.3276000000002, 
				10.4459000000003, 10.5814999999999, 10.7151999999998, 10.8589999999993, 
				10.9770999999993, 11.1487999999998, 11.2751999999993, 11.4159999999996, 
				11.5470999999998, 11.7013999999999, 11.8400000000001, 11.9778999999995, 
				12.1357000000004, 12.2686999999994, 12.3913999999993, 12.542400000001, 
				12.6807000000012, 12.8352000000003, 12.9618999999999, 13.1040000000005, 
				13.2717000000001, 13.3904999999999, 13.5542000000001, 13.6745999999996, 
				13.8351999999995), .Dim = c(100L, 2L), .Dimnames = list(c("10", 
						"20", "30", "40", "50", "60", "70", "80", "90", "100", "110", 
						"120", "130", "140", "150", "160", "170", "180", "190", "200", 
						"210", "220", "230", "240", "250", "260", "270", "280", "290", 
						"300", "310", "320", "330", "340", "350", "360", "370", "380", 
						"390", "400", "410", "420", "430", "440", "450", "460", "470", 
						"480", "490", "500", "510", "520", "530", "540", "550", "560", 
						"570", "580", "590", "600", "610", "620", "630", "640", "650", 
						"660", "670", "680", "690", "700", "710", "720", "730", "740", 
						"750", "760", "770", "780", "790", "800", "810", "820", "830", 
						"840", "850", "860", "870", "880", "890", "900", "910", "920", 
						"930", "940", "950", "960", "970", "980", "990", "1000"), c("parallel TRUE", 
						"parallel FALSE")))
plot(as.numeric(rownames(parstat)),parstat[,1],type='l',col="blue",ylim=c(0,14),xlab="delta (N in function arguments)")
lines(as.numeric(rownames(parstat)),parstat[,2],col="red")

