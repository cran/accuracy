\name{plot.accFirstLook}
\alias{plot.accFirstLook}
\alias{print.accFirstLook}
\title{ Generic methods for first look. }
\description{ Initial look objects are meant to printed or plotted }
\usage{
\method{plot}{accFirstLook}(x, single=FALSE, ask = dev.interactive(), ...)
\method{print}{accFirstLook}(x, ...)
}
\arguments{
  \item{x}{ perturb object to plotted or printed }
  \item{ask}{ ask to make plots}
  \item{\dots}{ additional parameters passed to xxx.default functions }
  \item{single}{ single plot only }
}

\value{
	see the related generic functions
}
\note{ 
 See \code{firstLook} for an example.}
\examples{
# see firstLook()
}
\author{
Micah Altman
\email{Micah\_Altman@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}
}


\seealso{ \code{\link{firstLook}}}

\keyword{misc}
\keyword{optimize}
\keyword{distribution}
