\name{R0}
\alias{R0}
\title{Calculates net reproduction, R0, according to a given set of rates Lx,fx and a fixed
 proportion female of births, \code{pfem}.}
\usage{
  R0(rates, pfem = 0.49)
}
\arguments{
  \item{rates}{rates needs to be a matrix with two columns,
  \code{Lx} and \code{Fx}. Here, \code{Lx} is the survival
  function integrated within each age interval and with a
  lifetable radix of 1. \code{Fx} is the fertility
  function, calculated as births/ person years of exposure.
  \code{sr} is the proportion female of births by each age
  of mother. \code{Fx} should simply contain zeros in ages
  with no fertility, OR, all vectors should be limited to
  reproductive ages. Place all vectors in a matrix, or
  specify rates as \code{cbind(Lx,Fx)}}

  \item{pfem}{the proportion female of births by each age
  of mother. Something like .49, .48, or (1/(2.05)).  This
  can either be specified as a single number, or it may be
  allowed to vary by age. For the later case, be sure to
  specify a value for each age (each row in \code{rates}).}
}
\value{
  the value of R0 for the given set of rates and proportion
  female of births.
}
\description{
  This function is only provided for the examples of
  \link{DecompContinuous}. It calculates the sum of the row
  products of \code{rates} multiplied by \code{pfem}.
}
\details{
  The main feature that functions need to have when
  specified for \code{DecompContinuous()} is that the rates
  must fit into a matrix, where presumably each column is a
  variable and each row is an age, as in \code{rates}.
  Really the decomposition function does not care how
  things are arranged in the matrix- the components of
  change matrix that is returned from
  \code{Decompcontinuous()} will be arranged in exactly the
  same way, so as long as you know how to interpret it, and
  your function can extract what it needs from the matrix,
  then it can be specified in any way. For this particular
  example function, \code{R0()}, \code{rates} must be
  specified with a variable in each column.
}
\note{
  It would also be possible to redefine the function to
  place \code{pfem} as an additional column in the rates
  matrix, which would allow this item to be decomposed too.
  Here it is specified separately in order to demonstrate
  passing on parameters to the function within
  \code{DecompContinuous()}.
}
\seealso{
  See Also as \code{\link{R0vec}}, the same function in the
  form necessary for \code{\link{DecompContinuousOrig}}
}

