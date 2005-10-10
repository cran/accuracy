\name{PTBdefault}
\alias{PTBdefault}
\alias{PTBdefaultfn}
\title{ Returns a default perturbation function for a given vector. }
\description{
 This function returns a function that can be used to  repeatedly peturb the vector given to it.
}
\usage{
PTBdefaultfn(vec, q = 0.99)
PTBdefault(vec, q=0.99)
}
\arguments{
  \item{vec}{ the vector which will be subject to perturbation}
  \item{q}{ for discrete vectors, the reclassification probability. For continuous vectors, perturbations will have 
	1-q relative uniform noise added. }
}

\value{
	A function that can be used to perturb the vector.
}
\details{
	For numeric discrete values, and ordered factors, \code{reclass.mat.diag} will
be used to generate the default reclassification matrix. For character vectors, unordered factors, and logicals,
\code{relass.mat.random} is used.
}

\examples{

x=1:100
rpx=replicate(100,x,simplify=FALSE) # perturb using the default method
# This produces equivalent results, but is faster, since reclass matrices are not recalculated on each
# replication
#fx=PTBdefault(x,q=.95) 
#rpx=replicate(100,fx(x),simplify=FALSE) 
matches = sapply(rpx,function(y)(sum(y==x)))	# how many matches to original vector
summary(matches)				# mean should be close to .95

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


\seealso{ \code{\link{perturb}},
\code{\link{PTBdiscrete}}, 
\code{\link{PTBus}}, 
\code{\link{reclass.mat.random}}}
\keyword{misc}
\keyword{optimize}
\keyword{distribution}
