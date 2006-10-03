\name{PTBi}
\alias{PTBi}
\alias{PTBms}
\alias{PTBmsb}
\alias{PTBmsbr}
\alias{PTBn}
\alias{PTBnc}
\alias{PTBns}
\alias{PTBnbr}
\alias{PTBnbrr}
\alias{PTBnsbr}
\alias{PTBnsbrr}
\alias{PTBu}
\alias{PTBuc}
\alias{PTBus}
\alias{PTBubr}
\alias{PTBubrr}
\alias{PTBusbr}
\alias{PTBusbrr}
\title{ functions to add random noise to a vector }
\description{
These functions add noise to a vector.
These are helper functions for \code{\link{perturb}}
They are used to select the type and magnitude of
noise applied to each vector in the data frame,
when running the perturbation sensitivity analysis. Use them only 
if they are substantively justified -- you can also supply custom
functions for use with this perturbation framework.
}

\usage{
	#
	# identity function
	#

	PTBi(x,size=1) 				# identify function, returns original vector
						# size is a dummy value, it is ignored
	#
	# minimal noise
	#

	PTBms(x, size=1)				# adds minimal, scaled, noise 
	PTBmsb(x, size=1, lbound=0, ubound=1)		# adds minimal, scaled, noise,
						# resulting vector truncated to bounds
	PTBmsbr(x, size=1, lbound=0, ubound=1)		# adds minimal, scaled  noise, out of-bounds elements
						# of vector resamples until within bounds

	#
	# Normally distributed noise
	#

	PTBn(x, size = 1)			# adds normally distributed noise
	PTBnc(x, size = 1)			# normally distributed noise, 
						# recentered to guarantee mean-zero disturbance
	PTBns(x, size = 1)			# adds normally distributed, scaled, noise

	PTBnbr(x, size = 1, lbound=0, ubound=1)# 
	PTBnsbr(x, size = 1, lbound=0, ubound=1)# adds normally (scaled)  distributed noise
						# out of-bounds elements
						# of vector resampled until within bounds

	PTBnbrr(x, size = 1, lbound=0, ubound=1)# 
	PTBnsbrr(x, size = 1, lbound=0, ubound=1)	# adds normally distributed,scaled noise
						# noise rescaled to bounds, and then 
						# out of bounds elements 
						# resampled until within bounds


	#
	# Uniformly distributed noise

	PTBu(x, size = 1)			# adds uniformly distributed noise
	PTBuc(x, size = 1)			# uniformly distributed noise, 
						# recentered to guarantee mean-zero disturbance
	PTBus(x, size = 1)			# adds uniformly distributed, scaled, noise

	PTBubr(x, size = 1, lbound=0, ubound=1)# 
	PTBusbr(x, size = 1, lbound=0, ubound=1)# adds uniformly, scaled,  distributed noise
						# out of-bounds elements
						# of vector resampled until within bounds

	PTBubrr(x, size = 1, lbound=0, ubound=1)# 
	PTBusbrr(x, size = 1, lbound=0, ubound=1)	# adds unirformly distributed,scaled noise
						# noise rescaled to bounds, and then 
						# out of bounds elements 
						# resampled until within bounds


}

\arguments{ 
  \item{x}{ vector to be perturbed }
  \item{size}{ magnitude of noise }
  \item{ubound}{ upper bound, can be a scalar or a vector of the same length as x }
  \item{lbound}{ lower bound, can be a scalar or a vector of the same length as x }
}

\details{

	Noise is randomly adds to each observation in the vector according to
	the chosen distribution function:
	- 'minimal' adds .Machine$double.ep) or subtracts .Machine$double.neg.eps
	- 'uniform' adds uniformly distributed random [-size/2,+size/2] 
	- 'normal' adds Normal(mean=0,stddev=size)
			 
	Noise is assymptotically mean zero. To guarantee mean-zero noise
	for a particular sample, use the "centered" helper functions, which
	adjust the noise vector to be approximately mean zero.

	The 'scaled' functions,  multiply each random deviate by $x_i$ for 
	each observation in x. This is useful if measurement error is
	heteroscedastic.

	Some helper functions can return a noisy vector guaranteed to be
	within given upper and lower bounds. This can be done through
	truncation, through resampling out-of-bounds observations, or
	by reducing the amount of noise added to entries close to the bounds.
	
}

\value{
	Returns a new vector x', with the same length as x.
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


\seealso{ \code{\link{perturb}}}
\examples{
	x=1:1000
	x.i = PTBi(x)
	x.m = PTBms(x)
	x.u = PTBu(x,size=.01)
	x.us = PTBus(x,size=.01)
	x.uc = PTBuc(x,size=.01)

	sum(x-x.i)		# should be 0, identity transform
	sum(which(x!=x.i))	# also 0

	sum(x-x.m)		# should be quite small, very close to 0
	length(which(x!=x.m))	# should be near 1000

	sum(x-x.u)		# should be small
	length(which(x!=x.m))	# should be near 1000

	sum(x-x.us)		# should be small
	length(which(x!=x.m))	# should be near 1000

	sum(x-x.uc)		# should be near 0
	length(which(x!=x.m))	# should be near 1000
	
}
\keyword{misc}
\keyword{optimize}
\keyword{distribution}
