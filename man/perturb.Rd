\name{sensitivty}
\alias{perturb}
\alias{sensitivity}
\title{ Data Perturbations based Sensitivity Analysis }
\description{
This function replicates an statistical analysis using slightly perturbed
versions of the original input data, and then analyzes the sensitivity of the
model to susch changes.  This can be used to draw attention to
        inferences that cannot be supported confidently given the current
        data, model, and algorithm/implementation,


}
\usage{
sensitivity(data, statistic, ..., ptb.R = 50, ptb.ran.gen = NULL , ptb.s = NULL, 
	summarize=FALSE)
perturb(...)
}
\arguments{
  \item{data}{ data matrix to be perturbed }
  \item{statistic}{ statistic model to be run on data }
  \item{\dots}{ additional arguments to statistical model}
  \item{ptb.R}{ number of replications for perturbation analysis }
  \item{ptb.ran.gen}{ a single function, or a vector of functions to be
		used to perturb each vector see \code{\link{PTBi}} }
  \item{ptb.s}{ a size, or vector of sizes, to be used in the vector perturbation functions}
  \item{summarize}{ if true, return a sensitivity summary, as would \code{summary(sensitivity())}. This 
reduces system memory use and can significantly speed up the analysis of large datasets. }

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
\note{  
If ptb.ran.gen is not specified,  then \code{PTBdefault} will be used, with q=ptb.s, for each
variable in the input data. 

Note "sensitivity" was  originally called "perturb". The name was changed to avoid a conflict
with another module, introduced afterwards.
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


\seealso{ See Also as \code{\link{PTBi}}, \code{\link{sensitivityZelig}}, \code{\link{PTBdefault}}, \code{\link{PTBdiscrete}}}
\examples{

# Examine the sensitivity of the GLM from Venables & Ripley (2002, p.189)
# as described in the glm module.
# 
# Perturb the two independent variables using +/- 0.025% 
#	(relative to the size of each observations)
# uniformly distributed noise. Dependent variable is not being modified.
# 
# Summary should show that estimated coefficients are not substantively affected by noise.

data(anorexia,package="MASS")
panorexia = sensitivity(anorexia, glm, Postwt ~ Prewt + Treat + offset(Prewt),
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
plongley = sensitivity(longley,lm,Employed~.) # defaults

# Alternatively, choose specific perturbation functions
#
#plongley = sensitivity(longley,lm,Employed~., ptb.R=100, 
#    ptb.ran.gen=c(PTBi, replicate(5,PTBus),PTBi), ptb.s=c(1,replicate(5,.001),1))

# summarizes range
sp=summary(plongley)
print(sp)
plot(sp) # plots rance of influence

# models with anova methods can also be summarized this way
anova(plongley) 

# plots different replications
plot(plongley) # plots the perturbed replications individually, pausing in between
plot(sp) # plots boxplots of the distribution of the coefficients under perturbation


# plots anova results (where applicable)
plot(anova(plongley))



# look in summary object to extract more ...
names(attributes(sp))

# print matrix of coefficients from all runs
coef= attr(sp,"coef.betas.m")
summary(coef)


# Example where model does not accept a dataset as an argument...

# MLE does not accept a dataset as an argument, so need to
# create a wrapper function
     library(stats4)
     mleD<-function(data,lld,...) {
           # construct LL function with embedded data
           f=formals(lld)
           f[1]=NULL
           ll <-function()  {
              cl=as.list(match.call())
              cl[1]=NULL
              cl$data=as.name("data")
              do.call(lld,cl)
           }
           formals(ll)=f

           # call mle
           mle(ll,...)
     }

     dat=as.data.frame(cbind( 0:10 , c(26, 17, 13, 12, 20, 5, 9, 8, 5, 4, 8) ))

     llD<-function(data, ymax=15, xhalf=6)
         -sum(stats::dpois(data[[2]], lambda=ymax/(1+data[[1]]/xhalf), log=TRUE))

     summary(sensitivity(dat, mleD,llD))



}
\keyword{misc}
\keyword{optimize}
