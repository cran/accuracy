\name{PTBdiscrete}
\alias{PTBdiscrete}
\title{Function to perturb vectors of discrete numeric values or factors, logicals, characters}
\description{
	This uses a reclassification matrix to perturb a vector discrete values.
}
\usage{
PTBdiscrete(x, r=NULL, q=0.99) # for factors, discrete numeric
}

\arguments{
  \item{x}{ the vector to be perturbed }
  \item{r}{ a reclassification matrix produced by  \code{reclass.mat.diag} or \code{reclass.mat.random}  or \code{reclassify} (from the separate perturb library) }
  \item{q}{ q-value for reclassification, if r is NULL}

}
\value{
	Returns a new vector x' with elements randomly reclassified according to r.
}

\details{ 
	This perturbs discrete vectors by reclassifying them at random based on a cumulative
probability reclassification matrix. 

	By default, if \code{r} is not supplied a default matrix will be supplied that has probability=\code{q} 
of keeping the same values, using the same rules as \code{PTBdefault}. 

As an alternative to classification for numeric vectors, see \code{PTBu} and related functions for methods that add random noise to a vector.
}

\references{
Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}
}

\author{
Micah Altman
\email{Micah\_Altman@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}
}

\seealso{ \code{\link{perturb}}, \code{\link{PTBn}}, \code{\link{PTBdefault}} ,
\code{\link{reclass.mat.random}}, \code{\link{reclass.mat.diag}} }

\examples{
	x=ceiling(runif(1:100)*3)		# vector taking on 3 discrete levels
	rx = reclass.mat.random(3,.95)	# reclassification matrix, prob of change  = .05
	rpx=replicate(100,PTBdiscrete(x,rx),simplify=FALSE) # 100 perturbations
	if (is.R()) {
     matches <-sapply(rpx,function(y)(sum(y==x)))	# how many matches to original vector
  } else {
     # Splus variation
     matches <-sapply(rpx,substitute(function(y)(sum(y==x))))
  }
	summary(matches)				# mean should be close to .95
}
\keyword{misc}
\keyword{optimize}
\keyword{distribution}


