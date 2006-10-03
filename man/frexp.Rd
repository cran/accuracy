\name{frexp}
\alias{frexp}
\title{Function to convert vector of floating-point numbers to fractional and integral components}
\description{
	The  frexp() function returns a matrix with the normalized fraction [.5,1) of the vector in the first
	column and the exponent (a power of two) in the second column. (If x is zero,
       then the normalized fraction is zero and zero is stored in the exponent.)
}
\usage{
	frexp(v)
}
\arguments{
  \item{v}{ vector of doubles }
}
\details{
  This is an R wrapper around the <math.h> frexp() function. 
}
\value{
Returns a 2xn matrix. Thethe normalized fraction [.5,1) of x is the first
	column and the exponent (a power of two) in the second column. (If x[i] is zero,
      then the normalized fraction of x[i] is zero and zero is stored in the exponent.)
}

\references{
Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}
}

\author{
Micah Altman
\email{Micah\_Altman\@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}
}


\examples{
x = runif(10)
y = frexp(x)
y
# this is now true by construction
x==y[,1] *2^y[,2]
\dontshow{
#unit test
if (is.R()){
  selfTest=get("frexpTest",envir=environment(frexp))
} else {
  selfTest=frexpTest
}
if (!selfTest()) {
  	warning("failed self-test")
}

}
}
\keyword{manip}
\keyword{misc}
\keyword{math}
