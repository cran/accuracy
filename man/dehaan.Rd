\name{dehaan}

\alias{dehaan}

\title{dehaan global optimality test}

\description{
Implements the de Haan test for identification of the global optimum of a likelihood
surface.
}


\usage{
   dehaan(llTest, llMax, pval=.05 )
}

\arguments{
  \item{llTest}{Vector of randomly generated likelihood values}
  \item{llMax}{Value of the likelihood function for a candidate global optimum}
  \item{pval}{p-value of the test}
}

\details{
\emph{dehaan} computes a (1-p) confidence interval for the global optimum of a
likelihood surface from a vector of user-supplied randomly chosen likelihood values.  
A user-supplied candidate maximum likelihood value is evaluated against confidence
interval.
}

\value{

\emph{dehaan=TRUE} if the candidate value is greater than the (1-p) confidence
interval for the true optimum.
}

\note{

The choice of the vector of user-supplied random likelihood values is somewhat of an
art.  The parameter space should be large enough to capture the true optimum, but not
so large as to include illegal or impractical parameter values.

The generated confidence interval depends on the number of random draws from the
parameter space.  The number should be large enough to allow the application of
asympotic theory.  A recommended number of random evaluations of the likelihood
function at a given set of parameter values is 500 or more, but depends on the size of
the parameter space.  Users are encouraged to experiment with more and less draws and
observe results.

}

\references{

Altman, M., J. Gill and M. P. McDonald.  2003.  \emph{Numerical Issues in Statistical
Computing for the Social Scientist}.  John Wiley \& Sons.
\url{http://www.hmdc.harvard.edu/numerical_issues/}

de Haan, L. 1981.  ``Estimation of the Minimum of a Function Using Order Statistics.''
\emph{Journal of the American Statistical Association} \bold{76}, 467-9.

Veall, M. R. 1990.  ``Testing for a Global Maximum in an Econometric Context.''
\emph{Econometrica} \bold{58} 1459-65.

}

\author{
Micah Altman
\email{Micah\_Altman@harvard.edu}
\url{http://www.hmdc.harvard.edu/micah_altman/}
, Michael P. Mcdonald
}

\examples{

  # The deHaan test is constructed as a maximum likelihood
  # test, with negative values for the likelihood.  The BOD problem
  # is a non-linear least squares minimization problem.  This test
  # is implemented using the negative of the sum of squares for consistency
  # with the deHaan framework of maximum likelihood.

  BOD <-
structure(list(Time = c(1, 2, 3, 4, 5, 7), demand = c(8.3, 10.3, 
19, 16, 15.6, 19.8)), .Names = c("Time", "demand"), row.names = c("1", 
"2", "3", "4", "5", "6"), class = "data.frame", reference = "A1.4, p. 270")
  stval<-expand.grid(A = seq(10, 100, 10), lrc = seq(.5, .8, .1))
  llfun<-function(A,lrc,BOD)
        -sum((BOD$demand - A*(1-exp(-exp(lrc)*BOD$Time)))^2)
  lls<-NULL 
  for (i in 1:nrow(stval))  {
     lls = rbind(lls, llfun(stval[i,1], stval[i,2],BOD))
  }
  fm1 <- nls(demand ~ A*(1-exp(-exp(lrc)*Time)),  
                data = BOD, start = c(A = 20, lrc = log(.35)))
  ss = -sum(resid(fm1)^2)
  dehaan(lls, ss)

\dontshow{
  if (is.R()){
    selfTest=get("dehaanSelftest",envir=environment(dehaan))
  } else {
    selfTest=dehaanSelftest
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

