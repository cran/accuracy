\name{perturb}
\alias{perturb}
\title{ Data Perturbations based Sensitivity Analysis }
\description{
This function replicates an statistical analysis using slightly perturbed
versions of the original input data, and then analyzes the sensitivity of the
model to susch changes.  This can be used to draw attention to
        inferences that cannot be supported confidently given the current
        data, model, and algorithm/implementation,


}
\usage{
perturb(data, statistic, ..., ptb.R = 1, ptb.ran.gen = PTBms, ptb.s = NULL)
}
\arguments{
  \item{data}{ data matrix to be perturbed }
  \item{statistic}{ statistic model to be run on data }
  \item{\dots}{ additional arguments to statistical model }
  \item{ptb.R}{ number of replications for perturbation analysis }
  \item{ptb.ran.gen}{ a single function, or a vector of functions to be
		used to perturb each vector see \code{\link{PTBi}} }
  \item{ptb.s}{ a size, or vector of sizes, to be used in the vector perturbation functions}
}
\details{

Sensitivity to numerical inaccuracy, and measurement error is very hard to 
measure formally. This empirical sensitivity tests draws 
attention to inferences that cannot be supported confidently given the current
data, model, and algorithm/implementation,

The empirical approach works by replicating
the original analysis while slightly perturbing the original  input data
in different ways. The sensitivity of the model estimates (e.g. 
estimated coefficients, standard errors and log-likelihoods) are then
summarized. 

The sensitivity analysis cannot be used to prove the accuracy of
a particular method, but is useful in drawing attention to
potential problems. Further experimentation and analysis may be
necessary to determine the specific cause of the problem:
numerical instability, statistical sensitivity to measurement error,
or ill-conditioning.

}
\value{
	Returns a list which contains the result  of each model run. Along
with attributes about the settings used in the perturbations.
	Use \code{summary} to summarize the results or extract a matrix of
 of the model parameters across the entire set of runs.
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


\seealso{ See Also as \code{\link{PTBi}}}
\examples{

# Examine the sensitivity of the GLM from Venables & Ripley (2002, p.189)
# as described in the glm module.
# 
# Perturb the two independent variables using +/- 0.025% 
#	(relative to the size of each observations)
# uniformly distributed noise. Dependent variable is not being modified.
# 
# Summary should show that estimated coefficients are not substantively affected by noise.

library(MASS)
data(anorexia)
panorexia = perturb(anorexia, glm, Postwt ~ Prewt + Treat + offset(Prewt),
     family=gaussian, 
    ptb.R=100, ptb.ran.gen=c(PTBi,PTBus,PTBus), ptb.s=c(1,.005,.005) )
summary(panorexia)

# Use classic longley dataset. The model is numerically unstable, 
# and much more sensitive to noise.  Smaller amounts of noise tremendously
# alter some of the estimated coefficients:
# 
# In this example we are not perturbing the dependent variable (employed) or 
# the year variable. So we assign then PTBi or NULL in ptb.ran.gen )

data(longley)
plongley = perturb(longley,lm,Employed~., ptb.R=100, 
    ptb.ran.gen=c(PTBi, replicate(5,PTBus),PTBi), ptb.s=c(1,replicate(5,.001),1))
sp=summary(plongley)

# summarizes range
print(sp)

# look in summary object to extract more ...
names(attributes(summary))

# print metrix of coefficients from all runs
coef= attr(sp,"coef.betas.m")
summary(coef)


}
\keyword{misc}
\keyword{optimize}
