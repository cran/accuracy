\name{trueRandom}
\alias{trueRandom}
\alias{resetSeed}
\alias{runifT}
\alias{initPool}
\title{ Function to return TRUE (not pseudo) random numbers, based on system and networked entropy collection. }
\description{
	This function makes use of hardware-generated entropy to supply 
true random numbers. Entropy collection is slow, so its primary use to
provide random numbers for cryptographic key generation and for setting 
the seed value for subsequent PRNG generation. 
}
\usage{
	runifT(n,min=0,max=1)
	trueRandom(n)
	resetSeed()
	initPool(size=512, hbok=TRUE, devrndok=TRUE) 
}
\arguments{
  \item{n}{ number of observations. If 'length(n) > 1', the length is taken to be the number required. }

 min,max: lower and upper limits of the distribution.
  \item{min,max}{ lower and upper limits of the distribution }
  \item{size}{ number of bytes of entropy to be retrieved into the local entropy pool}
  \item{hbok}{ use the "hotbits" entropy server }
  \item{devrndok}{ use /dev/random on Unix-derived systems }
}

\details{
These routines	provides true random numbers by accessing external entropy collectors.
Currently entropy can be retrieved from two different
sources. (1) The "Hotbits" server \url{http://www.fourmilab.ch/hotbits/} supplies
random bytes based on radioactive decay. (2) On Unix systems, the the  kernel gathers
environmental noise from device drivers and other sources into an entropy  pool, which
can be accessed through '/dev/random'.

By default, /dev/random will be preferred, if available, with hotbits used only
if /dev/random is not found. In addition, entropy will be requested from the
external sources in 512byte blocks (runifT will automatically request more
blocks if necessary to satisfy a request). Use initPool to change the source and
block size.

}

\value{
	Returns a vector of uniformly distributed true random numbers.
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
\warning{  This relies on external sources of true random numbers.
If a source of true random numbers is unavailable reset.seed() uses the clock value
and true.random() will return NULL.
	 }

\seealso{ See Also 
	\code{\link{runif}},
	\code{\link{.Random.seed}},
}
\examples{
if (.Platform$OS.type == "unix" || !is.null(nsl("www.cran.org"))) {
	# note, if you are using Windows, you must be on-line
	# to get to entropy generator
# reset the seed for runif() based on a true random value
resetSeed()
y=runif(100)

# generate true random values directly (may block for long periods if 
# if entropy pool is empty)

y=runifT(1000)
ty= table(trunc(5*y))
ty
chisq.test(ty)
} else {
	print ("Sorry, for entropy gathering on windows, requires a net connection.");
}

}
\keyword{misc}
\keyword{datagen}
\keyword{distribution}