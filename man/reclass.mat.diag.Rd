\name{reclass.mat.diag}
\alias{reclass.mat.diag}
\alias{reclass.mat.random}
\title{ Function to produce reclassification matrices}
\description{
This returns a cumulative probability matrix useful for reclassifying discrete variables with PTBdiscrete.
}
\usage{
reclass.mat.diag(n, q)
}
\arguments{
  \item{n}{ number of levels of factor, character vector, or discrete values}
  \item{q}{ probability of reclassifying to same sate }
}
\details{
  For reclass.mat.diag, the transition probability of changing from level i to level j is q if i==j, 0 if
  abs(i-j)>1 and q, 1-q/2 if abs(i-j)>1 and 1>i>n and 1-q otherwise

  For reclass.mat.random, the transition probability of changing from level i to level j is q if i==j, 
	 (1-q)/(n-1) otherwise.

}
\value{
	returns a matrix of cumulative probability distributions for reclassifying each level to
another level
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


\seealso{ \code{\link{perturb}}, \code{\link{reclass.mat.random}}}
\examples{
	x=ceiling(runif(1:100)*3)		# vector taking on 3 discrete levels
	rx = reclass.mat.random(3,.95)	# reclassification matrix, prob of change  = .05
	rpx=replicate(10,PTBdiscrete(x,rx),simplify=FALSE) # 100 perturbations
	matches = sapply(rpx,function(y)(sum(y==x)))	# how many matches to original vector
	summary(matches)				# mean should be close to .95
}
\keyword{misc}
\keyword{optimize}
\keyword{distribution}
