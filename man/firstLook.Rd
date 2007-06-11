\name{firstLook}
\alias{firstLook}
\title{Take a first look at a dataset}
\description{
	This takes a first look at the multivariate structure of a dataset
}
\usage{
	firstLook(x,forceType=TRUE)
}
\arguments{
  \item{x}{ a data frame}
  \item{forceType}{ whether to force non-numeric types for better exploration}
}
\details{
	This is for initial data exploration. Currently it uses PCA, parallel coordinates, and a heatmap of the correlation matrix.
}
\value{
	a firstlook object suitable for printing or plotting
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
	print(fl<-firstLook(swiss))
	plot(fl)
}
\keyword{manip}
\keyword{misc}
\keyword{math}
