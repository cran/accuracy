\name{plot.sensitivity}
\alias{plot.sensitivity}
\alias{summary.sensitivity}
\alias{print.sensitivity}
\alias{anova.sensitivity}
\alias{plot.sensitivity.summary}
\title{ Generic methods for perturbations. }
\description{
 Perturbations objects can be printed, summarized, and plotted .
}
\usage{
plot.sensitivity(x, ask=dev.interactive(), ...)
anova.sensitivity(object, ...)
print.sensitivity(x, quiet=TRUE, ...)
summary.sensitivity(object, ...)
}
\arguments{
  \item{x}{ perturb object to plotted or printed }
  \item{object}{ perturb object to be summarized, or presented as an anova}
  \item{\dots}{ additional parameters passed to xxx.default functions }
  \item{ask}{ whether to pause and ask between plots }
  \item{quiet}{ whether to print low-level details of perturb object}
}

\value{
	see the related generic functions
}
\note{ A perturb object can be summarized, printed, plotted, or summarized via anova if applicable. 
The idiom \code{plot(anova(perturb(...)))} will plot the anova summaries, as well. And the \code{plot(summary(perturb(...)))} will plot the sensitivity of the coefficients.
 See \code{perturb} for an example.}
\examples{
# see perturb()
}
\author{
Micah Altman
\email{Micah\_Altman@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}
}


\seealso{ \code{\link{perturb}}}

\keyword{misc}
\keyword{optimize}
\keyword{distribution}
