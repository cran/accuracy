\name{runifS}
\alias{trueRandom}
\alias{resetSeed}
\alias{runifT}
\alias{rnormS}
\alias{rstarMix}
\alias{runifS}
\title{ Function to return TRUE (not pseudo) random numbers, based on system and networked entropy collection. }
\description{
	This function makes use of hardware-generated entropy to supply 
true random numbers. Entropy collection is slow, so its primary use to
provide random numbers for cryptographic key generation and for setting 
the seed value for subsequent PRNG generation. 
}
\usage{
	runifT(n,min=0,max=1,maxTries=10,silent=TRUE)
	runifS(n,...,period=10000,maxTries=10,silent=TRUE)
	rnormS(n,...,period=10000,maxTries=10,silent=TRUE)
	rstarMix(n,...,period=10000,maxTries=10,silent=TRUE,rfunc="runif")
	resetSeed(maxTries=10,silent=TRUE,kind=NULL,normal.kind=NULL)
}
\arguments{
  \item{n}{ number of observations. If 'length(n) > 1', the length is taken to be the number required. }

  \item{min,max}{lower and upper limits of the distribution }
  \item{\dots}{pass on to runif}
  \item{period}{number of random numbers to generate before reseeding}
  \item{maxTries}{number of times to attempt to contact entropy source before giving yp}
  \item{silent}{whether to report failures of entropy source generators}
  \item{rfunc}{random deviate function to use}
  \item{kind}{sets the kind of random number generator as per \code{set.seed}}
  \item{normal.kind}{sets the kind of random number generator as per \code{set.seed}}
}

\details{
\code{runifS} returns uniform pseudo-random numbers, reseeding with a true random number every period draws.
\code{rnormS} returns normal pseudo-random numbers, reseeding with a true random number every period draws.
\code{rstarMix} returns random numbers from any specified random number deviant generating function, reseeding with a true random number every period draws. \code{rnormS} and \code{runifS} are simply wrappers around this function
\code{resetSeed} resets R's own PRNG seed, \code{set.seed}, supplying a true random number

\code{runifT} returns true random numbers in a uniform distribution. This can easily exhaust
entropy sources, so for larger samples, \code{runifS} is usually preferred.

These routines  provides true random numbers by accessing external entropy collectors.
Currently entropy can be retrieved from two different
sources. (1) The "Hotbits" server \url{http://www.fourmilab.ch/hotbits/} supplies
random bytes based on radioactive decay.
 (2) The bit server on random.org, and (3) On Unix systems, the the  kernel gathers
environmental noise from device drivers and other sources into an entropy  pool, which
can be accessed through '/dev/random'.

Entropy is retrieved in blocks from the sources (each source having a different preferred
block size). If the sources do not return a block within the R timeout value
(see \code{options("timeout")} entropy requests will be repeated up to maxTries times
per source. If no entropy is available, pseudo random numbers will be returned using \code{runif()}
(or a seed will be set based on \code{Sys.time}, for \code{resetSeed}, and a warning issued.

}

\value{
\code{runifS,runifT} 	Returns a vector of uniformly distributed true random numbers.
\code{resetSeed} returns the status of \code{set.seed()}

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

\seealso{ See Also 
	\code{\link{runif}},
	\code{\link{options}},
	\code{\link{.Random.seed}},
}
\examples{
# note, if you are using Windows, you must be on-line
# to get to entropy generator, will fall-back to pseudo-random
# numbers when off-line

\dontshow{
# only needed  to pass R CMD check on online windows
ow = options("warn"=0)
}

# reset the seed for runif() based on a true random value

resetSeed()
y=runif(100000)
ty= table(trunc(5*y))
print(ty)
if (is.R()) {
  chisq.test(ty)
} else {
  chisq.test(t(cbind(as.numeric(names(ty)),as.matrix(ty))))
}

\dontrun{
# generate true random values directly (may block for long periods if 
# if entropy pool is empty)
y=runifS(100000)
ty= table(trunc(5*y))
print(ty)
if (is.R()) {
  chisq.test(ty)
} else {
  chisq.test(t(cbind(as.numeric(names(ty)),as.matrix(ty))))
}

y=runifT(100)
ty= table(trunc(5*y))
print(ty);
chisq.test(ty)
}
\dontshow{
options(ow);
}

}
\keyword{misc}
\keyword{datagen}
\keyword{distribution}
