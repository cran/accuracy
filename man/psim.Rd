\name{psim}
\alias{psim}
\alias{setx.sensitivity}
\title{perform zelig simulations for perturbed models}
\description{
	This performs \code{sim()} for perturbed models.
}
\usage{
	psim(object,x=setx(object),...)
	 setx.sensitivity(obj,...)

}

\arguments{
  \item{object}{ an object of type perturb, as returned from \code{pzelig()}}
  \item{obj}{ an object of type perturb, as returned from \code{pzelig()}}
  \item{x}{ an object returned from \code{setx()} }
  \item{\dots}{ arguments to pass to default functions}
}

\value{
	PSim returns an object containing predictive simulations for each perturbation. 
	Setx returns a setx object. See \code{sim} and \code{setx} in the Zelig package.
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


\seealso{ See Also as \code{\link{perturb}},  \code{\link{PTBdefault}},  
	\code{\link[pkg:Zelig]{setx}},
	\code{\link[pkg:Zelig]{sim}},
		\code{\link{pzelig}} }
\examples{
	# see the example in pzelig()
}
\keyword{misc}
\keyword{optimize}
