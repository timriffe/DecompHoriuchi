
# Author: tim
###############################################################################

#' @title implementation of the decomposition algorithm of stepwise replacement  
#' @description This implements the algorithm described in Andreev et al (2002), with defaults set
#' to approximate their recommendations for replacement ordering and result averaging.
#' @details The \code{symmetrical} argument toggles whether or not we replace rates1 with rates2 (\code{FALSE}), 
#' or take the arithmetic average or replacement in both directions. \code{direction} refers to whether we go 
#' from the bottom up or top down, or take the arithmetic average of these when replacing vector elelements. 
#' Although the total difference will always sum correctly, the calculated contribution from individual components 
#' can vary greatly depending on the order in general. Defaults are set to symmetrically replace from the bottom 
#' up, per the authors' suggestion.
#' 
#' @param func A function specified by the user. This must be able to take the vectors \code{rates1} or 
#' \code{rates2} as its argument, and to return the value of the function, \code{y}, when evaluated for 
#' these rates. It may also have additional arguments, not to be decomposed.
#' @param rates1 vector of covariates to be passed on as arguments to \code{func()}. Covariates
#' can be in any order, as long as \code{func()} knows what to do with them. \code{rates1} is for time 1 
#' (or population 1).
#' @param rates2 is the same as \code{rates1} but for time/population 2.
#' @param symmetrical logical. default \code{TRUE} as recommended by authors. Shall we average the results of replacing 1 with 2 and 2 with 1?
#' @param direction. character. One of \code{"up"}, \code{"down"}, or \code{"both"}. Default \code{"up"}, as recommended by authors.
#' @param \dots optional parameters to pass on to \code{func()}.
#' 
#' @references E. Andreev, V. Shkolnikov, and A. Begun (2002) Algorithm for decomposition of differences between aggregate demographic measures and its application to life expectancies, healthy life expectancies, parity-progression ratios and total fertility rates. Demographic Research v7, n14
#' @export
#' @return a matrix of the variable effects that is organized in the same way as 
#' \code{rates1} and \code{rates2}.

stepwise_replacement <- function(func, rates1, rates2, symmetrical = TRUE, direction = "up",...){
	direction <- tolower(direction)
	stopifnot(direction %in% c("up","down","both"))
	
	up                   <- direction %in% c("up","both")
	down                 <- direction %in% c("down","both")
	N                    <- length(rates1)
	
	Rates1Mat            <- matrix(rates1, ncol = N + 1, nrow = N)
	Rates2Mat            <- matrix(rates2, ncol = N + 1, nrow = N)
	
	RM_1_2_up            <- matrix(ncol = N + 1, nrow = N)
	RM_1_2_down          <- RM_1_2_up
	RM_2_1_up            <- RM_1_2_up
	RM_2_1_down          <- RM_1_2_up
	
	# based on 1-> 2 upward
	r1ind                     <- lower.tri(Rates1Mat, TRUE)
	r2ind                     <- upper.tri(Rates1Mat)
	
	RM_1_2_up[r1ind]          <- Rates1Mat[r1ind]
	RM_1_2_up[r2ind]          <- Rates2Mat[r2ind]
	
	RM_1_2_down[r1ind[N:1, ]] <- Rates1Mat[r1ind[N:1, ]]
	RM_1_2_down[r2ind[N:1, ]] <- Rates2Mat[r2ind[N:1, ]]
	
	RM_2_1_up[r1ind]          <- Rates2Mat[r1ind]
	RM_2_1_up[r2ind]          <- Rates1Mat[r2ind]
	
	RM_2_1_down[r1ind[N:1, ]] <- Rates2Mat[r1ind[N:1, ]]
	RM_2_1_down[r2ind[N:1, ]] <- Rates1Mat[r2ind[N:1, ]]
	
	dec                       <- matrix(NA, nrow = N, ncol = 4)
	if (up){
		dec[, 1]              <- diff(apply(RM_1_2_up, 2, func, ...))
	}
	if (down){
		dec[, 2]              <- diff(apply(RM_1_2_down, 2, func, ...))
	}
	
	if (symmetrical){
		if (up){
			dec[, 3]          <- -diff(apply(RM_2_1_up, 2, func, ...))
		}
		if (down){
			dec[, 4]          <- -diff(apply(RM_2_1_down, 2, func, ...))
		}
	}
	
	
	dec_avg                   <- rowMeans(dec, na.rm = TRUE)
	dec_avg
}

#' @title an abridged lifetable based on M(x)
#' @description Implements the abridged lifetable formulas given in the supplementary material to Andreev et al (2012). An entire lifetable is calculated, but only lfie expectancy at birth is returned.
#' @details Chiang's a(x) is assumed in the following way: \eqn{a(0) = 0.07 + 1.7 * M(0)}, \eqn{a(1) = 1.6}, \eqn{a(\omega) = \frac{1}{M(\omega)}}, and all others are assumed at mid interval. The last age is assumed open. Everything else is pretty standard.
#' @param Mx numeric vector of abridged mortality rates.
#' @param Age integer, abrdged age lower bounds. 
#' @param radix numeric. Can be anything positive.
#' @references E. Andreev, V. Shkolnikov, and A. Begun (2002) Algorithm for decomposition of differences between aggregate demographic measures and its application to life expectancies, healthy life expectancies, parity-progression ratios and total fertility rates. Demographic Research v7, n14
#' 
#' E. Andreev, V. Shkolnikov, and A. Begun (2012) An Excel spreadsheet for the decomposition of a difference between two values of an aggregate demographic measure by stepwise replacement running from young to old ages. MPIDR Technical Report TR-2012-002
#' @export
#' @return numeric life expectancy at birth
LTabr <- function(Mx, Age = c(0,1,cumsum(rep(5,length(Mx)-2))),radix = 1e5){
	# based on lifetable formulas in Spreadsheet:
	# Andreev & Shkolnikov (2012) An Excel spreadsheet for the 
	# decomposition of a difference between two values of an 
	# aggregate demographic measure by stepwise replacement 
	# running from young to old ages 
	# MPIDR technical report 2012-002
	# attached file tr-2012-002-files.zip
	# containing spreadsheet
	# "Decomposition_replacement_from_young_to_old_ages(1).xls
	
	# TR: this is not HMD exactly, since it's an abridged lifetable. a(x) values
	# will all be different from HMD. Other potential differences are of lesser consequence
	# i.e. closeout at age 85, etc.
	N     <- length(Mx)
	w     <- c(diff(Age),5)
	
	ax    <- c(0.07 + 1.7 * Mx[1],.4,rep(.5,N - 3), 1/Mx[N])
	Nax   <- w * ax
	qx    <- (w * Mx) / (1 + w * (1 - ax) * Mx)
	qx[N] <- 1
	px    <- 1 - qx
	lx    <- c(radix, radix * (cumprod(px[-N])))
	dx    <- -diff(c(lx, 0))
	Lx    <- lx[-1] * w[-N] + dx[-N] * Nax[-N]
	Lx[N] <- ax[N] * lx[N]
	Tx    <- rev(cumsum(rev(Lx)))
	ex <- Tx / lx
	ex[1]
}

#' @title get life expectancy at birth from an (abridged)age-cause matrix
#' @description Given a matrix with abridged ages in rows and causes of death in columns, then calculate life expectancy at birth using \code{LTabr()}.
#' @details This assumes that the marginal row sums give all-cause mortality rates. Give an other category if you need to top-up to all-cause mortality. Do not include all-cause mortality itself!
#' @param Mxc numeric matrix
#' @return numeric life expectancy at birth
#' @export
Mxc2e0abr <- function(Mxc){
	Mx <- rowSums(Mxc)
	LTabr(Mx)
}
#' @title get life expectancy at birth from the vec of an age-cause matrix
#' @description Given a vector with abridged ages stacked within causes of death, assign its dimensions, take the age marginal sums using \code{Mxc2e0abr}, then calculate life expectancy at birth using \code{LTabr()}.
#' @details This assumes that the marginal row sums give all-cause mortality rates. Give an other category if you need to top-up to all-cause mortality. Do not include all-cause mortality itself! \code{length(Mxcvec)} must equal \code{prod(dim(Mxc))}. This function is meant to be fed to a generic decomposition function, such as \code{stepwise_replacement()}, or \code{DecompContinuousOrig()}.
#' @param Mxcvec numeric vector, \code{c(Mxc)}.
#' @param dims integer vector of length two, \code{c(nrow(Mxc),ncol(Mxc))}.
#' @param trans do we need to transpose in order to arrive back to an age-cause matrix?
#' @return numeric life expectancy at birth
#' @export
Mxc2e0abrvec <- function(Mxcvec, dims, trans = FALSE){
	dim(Mxcvec) <- dims
	if (trans){
		Mxcvec <- t(Mxcvec)
	}
	Mxc2e0abr(Mxcvec)
}

#' Year 2002 death rates by cause for US males in abridged age classes
#'
#' A matrix containing death rates for six causes (one of which is other) for abrdged age classes 0-85. Ages are labelled in rows, and causes in column names.
#'
#' @format A matrix with 19 rows and 6 columns

#' @source \url{https://www.demogr.mpg.de/en/projects_publications/publications_1904/mpidr_technical_reports/an_excel_spreadsheet_for_the_decomposition_of_a_difference_between_two_values_of_an_aggregate_4591.htm}
"Mxc1"

#' Year 2002 death rates by cause for England and Wales males in abridged age classes
#'
#' A matrix containing death rates for six causes (one of which is other) for abrdged age classes 0-85. Ages are labelled in rows, and causes in column names.
#'
#' @format A matrix with 19 rows and 6 columns

#' @source \url{https://www.demogr.mpg.de/en/projects_publications/publications_1904/mpidr_technical_reports/an_excel_spreadsheet_for_the_decomposition_of_a_difference_between_two_values_of_an_aggregate_4591.htm}
"Mxc2"

#' Comparison decomposition results by age and cause
#'
#' A matrix containing the contributions to the difference in life expectancy at birth between 2002 US males and England and Wales males. Ages (in rows) are in abridged categories, 0-85, and there are six causes, including other, in columns. The sum of the matrix is the difference in life expectancy at birth between the two populations. Values are based on symmetrical stepwise replacement from young to old ages only. This is just to make sure implementation is close.
#'
#' @format A matrix with 19 rows and 6 columns

#' @source \url{https://www.demogr.mpg.de/en/projects_publications/publications_1904/mpidr_technical_reports/an_excel_spreadsheet_for_the_decomposition_of_a_difference_between_two_values_of_an_aggregate_4591.htm}
"Compare"

