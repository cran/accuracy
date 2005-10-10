\name{pzelig}
\alias{pzelig}
\title{Perturbations-based Sensitivity Analysis of Zelig Models}
\description{
	This is a wrapper to perturb analyses launches through Zelig, with
	extension functions to summarize results.
}
\usage{
	pzelig(z, ptb.R = 50, ptb.ran.gen = NULL, ptb.s = NULL, explanatoryOnly=FALSE)
}
\arguments{
  \item{z} { an analysis output by \code{zelig} }
  \item{ptb.R}{ number of replications for perturbation analysis }
  \item{ptb.ran.gen}{ a single function, or a vector of functions to be
		used to perturb each vector see \code{\link{PTBi}} }
  \item{ptb.s}{ a size, or vector of sizes, to be used in the vector perturbation functions}
  \item{explanatoryOnly}{ a flag indicating whether explanatory variables should be perturbed by default}
}
\details{
	Uses \code{perturb} to perform sensistivity analysis on a \code{zelig} model.
}
\value{
	Returns a list which contains the result  of each model run. Along
with attributes about the settings used in the perturbations.
	Use \code{summary} to summarize the results or extract a matrix of
 of the model parameters across the entire set of runs.
}
\note{
	If ptb.ran.gen is not specified,  then \code{PTBdefault} will be used, with q=ptb.s, for each
	explanatory variable in the input data. By default, response variables will also 
	be perturbed, unless \code{explanatoryOnly} is true.
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


\seealso{ See Also as \code{\link{perturb}},  \code{\link{PTBdefault}},  \code{\link[pkg:Zelig]{zelig}},
		\code{\link{psim}} }
\examples{

# Sensitivity analysis of the classic longley data, as run through zelig.
#
# Note: Compare to the example documented under the main perturb page.

# example requires Zelig.

 if (library(Zelig,logical.return=TRUE)) {  
   data(longley)

   zelig.out=zelig(Employed~.,"ls",longley) # run longley regression
   # default settings
   perturb.zelig.out = pzelig(zelig.out)		   

   # Note: without zelig, the preceding would be:  perturb.out = perturb(longley,lm,Employed~.) 

   # plots the perturbed lm replications individually
   plot(perturb.zelig.out) 	  
   # summarizes the overall sensitivity of coefficient estimates
   summary(perturb.zelig.out)		  
   
   # since the "ls" model relies on "lm" which supports anova, can summarize
   # using anova
   anova(perturb.zelig.out)
   plot(anova(perturb.zelig.out))
   
   # Up to this point, we could have done this without Zelig... here's the fun part.
  
   # set values of explanatory variables
   setx.out=setx(perturb.zelig.out, Year=1955)       

   # Note: could also have used setx.out=setx(z,Year=1955)  instead of the line above,
   #   	   or  simply called psim() below without setx.out, which would default to using setx(perturb.zelig.out)

   # use simulation to predict value of explanatory variable based on
   # 
   sim.perturb.zelig.out = psim(perturb.zelig.out,setx.out)

   # this plots the profile of the predicted distribution of the explanatory variable, Employed
   # around the point wher Year=1955, and other explanatory variables are at their midpoints

   plot(sim.perturb.zelig.out)

   # this provides a summary of the predicted distribution of the explanatory variable

   summary(sim.perturb.zelig.out)
   
 }
}
\keyword{misc}
\keyword{optimize}