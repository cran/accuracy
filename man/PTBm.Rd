\name{PTBmu.gen}
\alias{PTBmu.gen}
\alias{PTBmn.gen}
\title{generator functions for multiple rounds of noise}
\description{
These functions generate a function that will apply multiple rounds of noise.
These are helper functions for \code{\link{perturb}}, and generate functions of the form of
\code{\link{PTBn}}.
They are used to select the type and magnitude of
noise applied to each vector in the data frame,
when running the perturbation sensitivity analysis. Use them only 
if they are substantively justified -- you can also supply custom
functions for use with this perturbation framework.
}

\usage{
	PTBmu.gen(reps=1)			# return a perturbation function that
	PTBmn.gen(reps=1)			# applies multiple rounds of noise
}

\arguments{ 
  \item{reps}{ number of rounds of noise to apply }
}

\details{
	Used to generate functions for use in \code{\link{perturb}}
}

\value{
	Returns a function of the form \code{\link{PTBn}} or 
	\code{\link{PTBi}}, with the modification that \emph{reps} 
	rounds of noise are applied for each perturbation.
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


\seealso{ \code{\link{perturb}}, \code{\link{PTBi}}}
\examples{
	x=1:1000
	f1=PTBmu.gen(); # should be roughly equivalent to PTBu() 
	x.u = f1(x,size=1)
	mean(x-x.u)  #should be small
	f2=PTBmu.gen(reps=100); # multiple disturbances tend to cancel eachother out
	x.u2 = f2(x,size=1)
	mean(x-x.u2)  #should be smaller 
}
\keyword{misc}
\keyword{optimize}
\keyword{distribution}
