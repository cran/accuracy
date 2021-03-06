\name{starr}
\alias{starr}

\title{staff global optimum test}

\description{ Implements the Starr test for identification of the global optimum of a
	likelihood surface. }

\usage{
 starr(betas, tol=.0001, dmethod="euclidean")
}

\arguments{
  \item{betas}{Vector of parameter values}
  \item{tol}{Tolerance distance between two parameter vectors to consider as "unique"
	optimum}
  \item{dmethod}{method used to compute distance between two parameter vectors}
}

\details{
	\emph{starr} computes the probability that a local optimum, which may or may
not be the global optimum, of a likelihood function has not been observed.  The
probability is generated by counting the observed number of ``basins of attraction'' -
starting values that lead to an local optimum.
}

\value{ 
\emph{starr} is a class "double" with value equal to the probability than a
	local optimum is unobserved.
}

\note{
\emph{starr} is given a vector of \emph{betas} parameter values identifying local optima resultant of a search algorithm, such as \emph{mle} or \emph{nls}.  The starting values used to generate \emph{betas} can be a grid, or for larger n-dimenional parameter spaces, may be randomly chosen.  The parameter values identifying a local optimum are passed to \emph{starr} for each set of starting values that lead to that local optimum.
}



\references{
Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}

Finch, S. J., Mendell, N. R., Thode, H. C., Jr. 1989. ``Probabilistic Measures of
Adequacy of a Numerical Search for a Global Maximum.'' \emph{Journal of the American
Statistical Association} \bold{84}, 1020-3.

Starr, N. 1979. ``Linear Estimation of the Probability of Discovering a New Species.''
\emph{Annals of Statistics} \bold{7}, 644-52.
}

\author{
Micah Altman
\email{Micah\_Altman@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}
, Michael McDonald
}


\examples{

	x=rbind(c(1,1,1), c(1,2,1), c(1,1.1,1), c(1,2,1), c(3,4,5));
	starr(rbind(1,1,2,2));

      #BODmodel
        
      BOD <-
structure(list(Time = c(1, 2, 3, 4, 5, 7), demand = c(8.3, 10.3, 
19, 16, 15.6, 19.8)), .Names = c("Time", "demand"), row.names = c("1", 
"2", "3", "4", "5", "6"), class = "data.frame", reference = "A1.4, p. 270")

       stval=expand.grid(A = seq(10, 100, 10), lrc = seq(.5, .8, .1))
       llfun<-function(A,lrc,BOD)
             -sum((BOD$demand - A*(1-exp(-exp(lrc)*BOD$Time)))^2)
       lls=NULL 
       for (i in 1:nrow(stval))  {
          lls = rbind(lls, llfun(stval[i,1], stval[i,2],BOD))
       }
       fm1 <- nls(demand ~ A*(1-exp(-exp(lrc)*Time)),  
                     data = BOD, start = c(A = 20, lrc = log(.35)))
       ss = -sum(resid(fm1)^2)
       dehaan(lls, ss)
       
       llb=NULL
       for (i in 1:nrow(stval))  {
	llb = rbind(llb,coef( 
  	nls(demand ~ A*(1-exp(-exp(lrc)*Time)),  
                     data = BOD, start = c(A=stval[i,1], lrc = stval[i,2]))))
      }
       starr(llb)
       
\dontshow{


if (is.R()){
  selfTest=get("starrSelftest",envir=environment(starr))
} else {
  selfTest=starrSelftest
}
if (!selfTest()) {
  	warning("failed self-test")
}
}

}
\keyword{manip}
\keyword{misc}
\keyword{math}
\keyword{htest}
