\name{HTML.perturbS}
\alias{HTML.perturbS}
\alias{HTML.perturb.sim.summary}
\alias{HTML.perturbAnova}
\title{HTML methods for perturbation analysis}
\description{
	HTML pretty printing methods for use rith R2HTML package
}
\usage{
	HTML.perturbS(x,...)
	HTML.perturb.sim.summary(x,...)
	HTML.perturbAnova(x,...)

}

\arguments{
  \item{x}{ a perturb summary object} 
  \item{\dots}{ arguments to pass to default functions}
}

\value{
	None
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


\seealso{ See Also as \code{\link{perturb}},  \code{\link[pkg:R2HTML]{HTML}}}
\examples{
    \dontrun{
library(R2HTML)
sp=summary(perturb(longley,lm,Employed~.))
HTML(sp,file="") # HTML Output
    }
}
\keyword{misc}
\keyword{optimize}
